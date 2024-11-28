/*
gcc cov_mat_an.c -lmseed -lfftw3 -llapack -lrefblas -lgfortran -lm
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <glob.h>

#include <libmseed.h>
#include <fftw3.h>

typedef struct
{
        double real;
        double imag;
} complex;

int sub_len=48;
int nsub_1avg=100;
int sampling_rate_0=200;
int sampling_rate=20;
double freq_min=0.01;
double freq_max=10;
int nsta_min=3;

int sec_1d=86400;

double delta;
int factor,npts_1d_0,npts_2d_0;
int sub_npts,sub_shift_npts;
int navg_1d,avg_shift_npts_0,avg_npts_0,avg_npts;

double dtrnd_denom,*dtrnd_dvtn;

int nfreq;
double *fftw_input;
fftw_complex *fftw_output;
fftw_plan plan;
int fmin_idx,fmax_idx,nfreq_band;

struct tm cur_day_bdt,next_day_bdt;
time_t cur_day_sct,next_day_sct;
int nday,navg;

int nsta=0;
char (*stnm)[5];
double *data;
int *mask;
int *day_gap;
int *sta_avg;
double *data_dtrnd;
double *data_avg;
complex *sta_vctr;
complex *ap;
complex *cov_mat;

double *w,*rwork;
complex *z,*work;
int *iwork;
int lwork=-1,lrwork=-1,liwork=-1;

int nsta_avg,npack_avg;

int glob_flags=0;
int8_t splitversion=0;
uint32_t mstl_flags=0;
int8_t verbose=0;
int8_t truncate=0;
int8_t freeprvtptr=1;

void zhpr_(char *uplo,int *n,double *alpha,complex *x,int *incx,complex *ap);
void zhpevd_(char *jobz,char *uplo,int *n,complex *ap,double *w,complex *z,int *ldz,complex *work,int *lwork,double *rwork,int *lrwork,int *iwork,int *liwork,int *info);

void init_global(void)
{
	int sub_shift,sub_ovrlp;
	int avg_len,nsub_1shift_avg,avg_shift;

	int i;
	
	delta=(double)1/sampling_rate;
	factor=sampling_rate_0/sampling_rate;
	npts_1d_0=sec_1d*sampling_rate_0;
	npts_2d_0=2*npts_1d_0;

	sub_npts=sub_len*sampling_rate;
	sub_shift=sub_len/2;
	sub_ovrlp=sub_len-sub_shift;
	sub_shift_npts=sub_shift*sampling_rate;

	avg_len=sub_shift*nsub_1avg;
	nsub_1shift_avg=nsub_1avg/2;
	avg_shift=sub_shift*nsub_1shift_avg;
	navg_1d=sec_1d/avg_shift;
	avg_shift_npts_0=avg_shift*sampling_rate_0;
	avg_npts_0=(avg_len+sub_ovrlp)*sampling_rate_0;
	avg_npts=(avg_len+sub_ovrlp)*sampling_rate;

	dtrnd_denom=(double)(avg_npts_0-1)*avg_npts_0*(avg_npts_0+1)/12;
	dtrnd_dvtn=(double *)malloc(avg_npts_0*sizeof(double));
	for(i=0;i<avg_npts_0/2;i++)
	{
		dtrnd_dvtn[avg_npts_0/2+i]=(double)(2*i+1)/2;
		dtrnd_dvtn[(avg_npts_0/2-1)-i]=-dtrnd_dvtn[avg_npts_0/2+i];
	}

	nfreq=sub_npts/2+1;
	fftw_input=fftw_alloc_real(sub_npts);
	fftw_output=fftw_alloc_complex(nfreq);
	plan=fftw_plan_dft_r2c_1d(sub_npts,fftw_input,fftw_output,FFTW_MEASURE);
	fmin_idx=10;
	fmax_idx=11;
	nfreq_band=(fmax_idx+1)-fmin_idx;

	mstl_flags|=MSF_UNPACKDATA;
}

void term_global(void)
{
	free(dtrnd_dvtn);
	
	fftw_free(fftw_input);
	fftw_free(fftw_output);
	fftw_destroy_plan(plan);

	free(stnm);
	free(data);
	free(mask);
	free(sta_avg);
	free(data_dtrnd);
	free(data_avg);
	free(sta_vctr);
	free(ap);
	free(cov_mat);

	free(w);
	free(z);
	free(work);
	free(rwork);
	free(iwork);
}

void init_date(char *begin_str,char *end_str)
{
	struct tm begin_bdt,end_bdt; //broken-down time
	time_t begin_sct,end_sct; //simple calendar time

	//begin date
	strptime(begin_str,"%Y-%m-%d",&begin_bdt);
	begin_bdt.tm_hour=0;
	begin_bdt.tm_min=0;
	begin_bdt.tm_sec=0;
	cur_day_bdt=begin_bdt;

	//end date
	strptime(end_str,"%Y-%m-%d",&end_bdt);
	end_bdt.tm_hour=0;
	end_bdt.tm_min=0;
	end_bdt.tm_sec=0;

	//no. of total days
	begin_sct=timegm(&begin_bdt);
	cur_day_sct=begin_sct;
	end_sct=timegm(&end_bdt);
	nday=(end_sct-begin_sct)/sec_1d+1;
	navg=nday*navg_1d;

	//next day
	next_day_sct=cur_day_sct+sec_1d;
	gmtime_r(&next_day_sct,&next_day_bdt);
}

void init_stnm(void)
{
	char stnm_file[9]="sta_meta";
	FILE *fp;
	char dummy_stnm[5];

	int i;

	fp=fopen(stnm_file,"r");

	//1. just read all to count station number
	while(fscanf(fp,"%s\n",dummy_stnm)!=EOF)
	{
		nsta++;
	}
	stnm=(char (*)[5])malloc(nsta*5*sizeof(char));

	//2. read again to store station names
	fseek(fp,0,SEEK_SET);
	for(i=0;i<nsta;i++)
	{
		fscanf(fp,"%s\n",stnm[i]);
	}

	data=(double *)malloc(nsta*npts_2d_0*sizeof(double)); //nsta x npts_2d_0
	mask=(int *)malloc(nsta*npts_2d_0*sizeof(int));
	day_gap=(int *)malloc(2*nsta*sizeof(int));
	sta_avg=(int *)malloc(nsta*sizeof(int));
	data_dtrnd=(double *)malloc(nsta*avg_npts_0*sizeof(double)); //nsta x avg_npts_0
	data_avg=(double *)malloc(nsta*avg_npts*sizeof(double)); //nsta x avg_npts
	sta_vctr=(complex *)malloc(nsta*nfreq_band*sizeof(complex)); //nsta x nfreq_band
	ap=(complex *)malloc(nsta*(nsta+1)/2*sizeof(complex));
	cov_mat=(complex *)malloc(nfreq_band*nsta*(nsta+1)/2*sizeof(complex)); //nfreq_band x nsta*(nsta+1)/2

	fclose(fp);
}

void init_eig(void)
{
	char jobz='V',uplo='U';
	int n=nsta;
	complex *dummy_ap=malloc(nsta*(nsta+1)/2*sizeof(complex));
	int ldz=n,info;

	w=malloc(nsta*sizeof(double));
	z=malloc(ldz*n*sizeof(complex));
	work=malloc(1*sizeof(complex));
	rwork=malloc(1*sizeof(double));
	iwork=malloc(1*sizeof(int));

	zhpevd_(&jobz,&uplo,&n,dummy_ap,w,z,&ldz,work,&lwork,rwork,&lrwork,iwork,&liwork,&info);

	lwork=work[0].real;
	work=malloc(lwork*sizeof(complex));
	lrwork=rwork[0];
	rwork=malloc(lrwork*sizeof(double));
	liwork=iwork[0];
	iwork=malloc(liwork*sizeof(int));

	free(dummy_ap);
}

int merge_seg(MS3TraceList *mstl,time_t day_sct,double *data,int *mask)
{
	MS3TraceID *tid;
	MS3TraceSeg *seg;
	nstime_t prev_endtime=0;
	int data_loff,seg_loff;
	int data_roff,seg_roff;
	int npts;

	int i;

	tid=mstl->traces.next[0];
	while(tid)
	{
		seg=tid->first;
		while(seg)
		{
			if(seg->endtime<prev_endtime)
			{
				//contained segment: just pass
				seg=seg->next;
				continue;
			}

			//convert segment integer to double
			if(mstl3_convertsamples(seg,'d',truncate))
			{
				fprintf(stderr,"Error in coverting integer segment to double\n");
				return 1;
			}

			data_loff=round((seg->starttime/1E9-day_sct)*sampling_rate_0);
			seg_loff=0;
			if(data_loff<0)
			{
				//segment left part lies outside data
				seg_loff=-data_loff;
				data_loff=0;
			}

			data_roff=round(((day_sct+sec_1d)-seg->endtime/1E9)*sampling_rate_0)-1;
			seg_roff=0;
			if(data_roff<0)
			{
				//segment right part lies outside data
				seg_roff=-data_roff;
				data_roff=0;
			}

			npts=seg->numsamples-seg_loff-seg_roff;

			for(i=0;i<npts;i++)
			{
				if(mask[data_loff+i]==0)
				{
					//data is empty -> fill it with segment
					data[data_loff+i]=*((double *)seg->datasamples+seg_loff+i);
					mask[data_loff+i]=1;
				}

				else if(mask[data_loff+i]==1)
				{
					//data already occpuies -> check mismatch between segment and data
					if(data[data_loff+i]!=*((double *)seg->datasamples+seg_loff+i))
					{
						//mismatch between segment and data
						mask[data_loff+i]=2;
					}
				}
			}

			prev_endtime=seg->endtime;
			seg=seg->next;
		}

		tid=tid->next[0];
	}

	return 0;
}

int read_mseed(char *stnm,time_t day_sct,double *data,int *mask)
{
	struct tm day_bdt;
	char pattern[68];
	glob_t matches;
	MS3TraceList *mstl=NULL;

	int status;
	int i=0;

	gmtime_r(&day_sct,&day_bdt);

	//construct pattern
	sprintf(pattern,"/home/tllc46/48NAS2/symbolic.jeju/05.%s/HHZ/05.%s.HHZ.%d.%03d*",stnm,stnm,1900+day_bdt.tm_year,1+day_bdt.tm_yday);

	//glob pattern match
	status=glob(pattern,glob_flags,NULL,&matches);
	if(status==GLOB_NOMATCH)
	{
		fprintf(stderr,"No matches at all\n");
		return 1;
	}
	else if(status)
	{
		fprintf(stderr,"Error in glob pattern matching\n");
		return -1;
	}

	//read mseed file to trace list
	while(i<matches.gl_pathc)
	{
		if(ms3_readtracelist(&mstl,matches.gl_pathv[i],NULL,splitversion,mstl_flags,verbose))
		{
			fprintf(stderr,"Error in reading mseed file to trace list\n");
			return -1;
		}
		i++;
	}

	//initialize mask array
	memset(mask,0,npts_1d_0*sizeof(int));

	//merge segments
	if(0<merge_seg(mstl,day_sct,data,mask))
	{
		fprintf(stderr,"Error in merging segments\n");
		return -1;
	}

	mstl3_free(&mstl,freeprvtptr);

	return 0;
}

int read_mseeds(time_t day_sct,double *data,int *mask,int *day_gap)
{
	int i;

	for(i=0;i<nsta;i++)
	{
		if(day_gap[i]=(read_mseed(stnm[i],day_sct,data+i*npts_2d_0,mask+i*npts_2d_0))==-1)
		{
			fprintf(stderr,"Error in reading %s data\n",stnm[i]);
			return 1;
		}
	}

	return 0;
}

void detrend(double *data)
{
	int sta_idx;
	double mean;
	double dtrnd_numer;
	double dtrnd_slope;

	int i,j;

	for(i=0;i<nsta_avg;i++)
	{
		sta_idx=sta_avg[i];
		mean=0;
		dtrnd_numer=0;

		for(j=0;j<avg_npts_0;j++)
		{
			mean+=data[sta_idx*npts_2d_0+j];
		}
		mean=mean/avg_npts_0;

		for(j=0;j<avg_npts_0;j++)
		{
			dtrnd_numer+=dtrnd_dvtn[j]*(data[sta_idx*npts_2d_0+j]-mean);
		}

		dtrnd_slope=dtrnd_numer/dtrnd_denom;

		for(j=0;j<avg_npts_0;j++)
		{
			data_dtrnd[i*avg_npts_0+j]=data[sta_idx*npts_2d_0+j]-mean-dtrnd_slope*dtrnd_dvtn[j];
		}
	}
}

void decimate(void)
{
	//just desampling... different from SAC's decimate and same as obspy's decimate with no_filter=True
	int i,j;

	for(i=0;i<nsta_avg;i++)
	{
		for(j=0;j<avg_npts;j++)
		{
			data_avg[i*avg_npts+j]=data_dtrnd[i*avg_npts_0+j*factor];
		}
	}
}

void cov_mat_avg(void)
{
	char uplo='U';
	int n=nsta_avg;
	double alpha=1;
	int incx=nfreq_band;

	int i,j,k;

	npack_avg=nsta_avg*(nsta_avg+1)/2;
	memset(cov_mat,0,nfreq_band*npack_avg*sizeof(complex));
	for(i=0;i<nsub_1avg;i++) //for each sub window
	{
		for(j=0;j<nsta_avg;j++) //for each station
		{
			for(k=0;k<sub_npts;k++)
			{
				fftw_input[k]=data_avg[i*sub_shift_npts+j*avg_npts+k];
			}

			fftw_execute(plan);

			for(k=0;k<nfreq_band;k++)
			{
				sta_vctr[j*nfreq_band+k].real=fftw_output[k+fmin_idx][0];
				sta_vctr[j*nfreq_band+k].imag=fftw_output[k+fmin_idx][1];
			}
		}

		for(j=0;j<nfreq_band;j++)
		{
			memset(ap,0,npack_avg*sizeof(complex)); //because zhpr_ function cumulates ap
			zhpr_(&uplo,&n,&alpha,sta_vctr+j,&incx,ap);

			for(k=0;k<npack_avg;k++)
			{
				cov_mat[j*npack_avg+k].real+=ap[k].real;
				cov_mat[j*npack_avg+k].imag+=ap[k].imag;
			}
		}
	}

	for(i=0;i<nfreq_band*npack_avg;i++)
	{
		cov_mat[i].real=cov_mat[i].real/nsub_1avg;
		cov_mat[i].imag=cov_mat[i].imag/nsub_1avg;
	}
}

void calculate_spcwdth(void)
{	
	int gap_flag;
	char jobz='V',uplo='U';
	int ldz,info;

	int i,j,k;

	//for(i=0;i<navg_1d;i++) //for each average window
	for(i=0;i<2;i++)
	{
		nsta_avg=0;
		for(j=0;j<nsta;j++) //for each station
		{
			gap_flag=0;
			for(k=0;k<avg_npts_0;k++)
			{
				if(mask[j*npts_2d_0+i*avg_shift_npts_0+k]!=1)
				{
					//average window contains gap
					gap_flag=1;
					break;
				}
			}

			if(!gap_flag)
			{
				sta_avg[nsta_avg]=j;
				nsta_avg++;
			}
		}

		if(nsta_avg<nsta_min)
		{
			fprintf(stderr,"%02d/%d average window: traces are too few\n",i+1,navg_1d);
			continue;
		}

		detrend(data+i*avg_shift_npts_0);
		decimate();
		cov_mat_avg();

		ldz=nsta_avg;
		for(j=0;j<nfreq_band;j++)
		{
			zhpevd_(&jobz,&uplo,&nsta_avg,cov_mat+j*npack_avg,w,z,&ldz,work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
		}
	}
}

void next_day(void)
{
	int i;
	
	//shift day
	cur_day_sct=next_day_sct;
	next_day_sct+=sec_1d;
	gmtime_r(&cur_day_sct,&cur_day_bdt);
	gmtime_r(&next_day_sct,&next_day_bdt);

	for(i=0;i<nsta;i++)
	{
		if(!day_gap[nsta+i])
		{
			//shift data and mask array
			memmove(data+i*npts_2d_0,data+i*npts_2d_0+npts_1d_0,npts_1d_0*sizeof(double));
			memmove(mask+i*npts_2d_0,mask+i*npts_2d_0+npts_1d_0,npts_1d_0*sizeof(int));
		}
	}

	memmove(day_gap,day_gap+nsta,nsta*sizeof(int));
}

int main(int argc,char **argv)
{
	char info[20];
	int nsta_day;

	int i,j,status;
	
	init_global();
	init_date("2014-04-10","2014-04-10");
	init_stnm();
	init_eig();

	if(read_mseeds(cur_day_sct,data,mask,day_gap))
	{
		fprintf(stderr,"Error in reading data\n");
		return 1;
	}

	//main loop
	for(i=0;i<nday;i++)
	{
		strftime(info,20,"%Y-%m-%d starting",&cur_day_bdt);
		printf("%s\n",info);

		//read next day mseed file to trace list
		if(read_mseeds(next_day_sct,data+npts_1d_0,mask+npts_1d_0,day_gap+nsta))
		{
			fprintf(stderr,"Error in reading data\n");
			return 1;
		}

		nsta_day=0;
		for(j=0;j<nsta;j++)
		{
			if(!day_gap[j])
			{
				nsta_day++;
			}
		}

		if(nsta_day<nsta_min)
		{
			fprintf(stderr,"traces are too few\n");
			/*for(j=0;j<navg_1d;j++)
			{
				fwrite(psd_gap,sizeof(double),nfreq_bin,fp);
			}*/
			next_day();
			continue;
		}

		calculate_spcwdth();

		next_day();
	}

	term_global();
	return 0;
}
