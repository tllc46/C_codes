#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <glob.h>

#include <libmseed.h>
#include <fftw3.h>

int sec_1d=86400;
int sampling_rate=200;
int ppsd_length=3600;
double overlap=0.5;
int period_smoothing_width_octaves=1;
double period_step_octaves=0.125;

int npts_1d;
double delta;
int nstep,nseg_1d,nfft,nhop,nwin,nfreq;
int nfreq_no_zero,nfreq_bin;
int *right_idx,*left_idx;
double *deviation;
int scaling_factor=2;
int nday,nseg;

int glob_flags=0;
int8_t splitversion=0;
uint32_t mstl_flags=0;
int8_t verbose=0;
int8_t truncate=0;
int8_t freeprvtptr=1;

void init_global(void)
{
	int seg_len,step,win_len;
	int exponent;
	double *frequency_log2;
	double right_exp,left_exp;

	int i,jr=0,jl=0;

	npts_1d=sec_1d*sampling_rate;
	delta=(double)1/sampling_rate;
	seg_len=sampling_rate*ppsd_length;
	step=ppsd_length*(1-overlap);
	nseg_1d=sec_1d/step;
	nstep=step*sampling_rate;

	//13 segments overlapping by 75%
	//1 segment + 25% * 12 segments
	win_len=ppsd_length*sampling_rate/4;

	exponent=log2(win_len);
	nfft=exp2(exponent); //previous power of 2
	nhop=nfft/4; //75% overlap -> 25% hop
	nwin=(seg_len-nfft)/nhop+1;
	nfreq=nfft/2+1;

	//log2 frequency
	nfreq_no_zero=nfreq-1;
	frequency_log2=(double *)malloc(nfreq_no_zero*sizeof(double));
	for(i=0;i<nfreq_no_zero;i++)
	{
		frequency_log2[i]=log2(i+1);
	}

	//frequency bin and smoothing window index
	nfreq_bin=(exponent-1)/period_step_octaves+1;
	right_idx=(int *)malloc(nfreq_bin*sizeof(int));
	left_idx=(int *)malloc(nfreq_bin*sizeof(int));
	for(i=0;i<nfreq_bin;i++)
	{
		right_exp=i*period_step_octaves+0.5*period_smoothing_width_octaves;
		while(frequency_log2[jr]<=right_exp && jr<nfreq_no_zero)
		{
			jr++;
		}
		right_idx[i]=jr;

		left_exp=i*period_step_octaves-0.5*period_smoothing_width_octaves;
		while(frequency_log2[jl]<left_exp && jl<nfreq_no_zero)
		{
			jl++;
		}
		left_idx[i]=jl;
	}

	//initialize deviation
	deviation=(double *)malloc(nfft*sizeof(double));
	for(i=0;i<nfft/2;i++)
	{
		deviation[i+nfft/2]=0.5*(2*i+1)*delta;
		deviation[-i+(nfft/2-1)]=-deviation[i+nfft/2];
	}

	mstl_flags|=MSF_UNPACKDATA;

	free(frequency_log2);
}

void term_global(void)
{
	free(right_idx);
	free(left_idx);
	free(deviation);
}

void init_date(char *begin_str,char *end_str,time_t *cur_day_smpl,struct tm *cur_day_tm)
{
	struct tm begin_tm,end_tm;
	time_t begin_smpl,end_smpl;

	strptime(begin_str,"%Y-%m-%d",&begin_tm);
	begin_tm.tm_hour=0;
	begin_tm.tm_min=0;
	begin_tm.tm_sec=0;
	*cur_day_tm=begin_tm;

	strptime(end_str,"%Y-%m-%d",&end_tm);
	end_tm.tm_hour=0;
	end_tm.tm_min=0;
	end_tm.tm_sec=0;

	begin_smpl=timegm(&begin_tm);
	*cur_day_smpl=begin_smpl;
	end_smpl=timegm(&end_tm);
	nday=(end_smpl-begin_smpl)/sec_1d+1;
	nseg=nday*nseg_1d;
}

void init_npy(FILE* fp)
{
	char magic[7]="\x93" "NUMPY";
	char major=1,minor=0;
	unsigned short header_len=118;
	char *header=(char *)malloc((header_len+1)*sizeof(char));
	double array[3];
	int status;

	memset(header,' ',header_len);

	status=sprintf(header,"{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d), }",nseg,nfreq_no_zero);

	header[status]=' ';
	header[header_len-1]='\n';
	header[header_len]='\0';

	fputs(magic,fp);

	fwrite(&major,sizeof(char),1,fp);
	fwrite(&minor,sizeof(char),1,fp);
	fwrite(&header_len,sizeof(unsigned short),1,fp);

	fputs(header,fp);
}

int merge_seg(MS3TraceList *mstl,time_t day_smpl,double *data,int *mask)
{
	MS3TraceID *tid;
	MS3TraceSeg *seg;
	nstime_t prev_endtime=0;
	int b_offset_data,b_offset_seg;
	int e_offset_data,e_offset_seg;
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
				prev_endtime=seg->endtime;
				seg=seg->next;
				continue;
			}

			//convert segment integer to double
			if(mstl3_convertsamples(seg,'d',truncate))
			{
				fprintf(stderr,"Error in coverting integer segment to double\n");
				return 1;
			}

			b_offset_data=(seg->starttime/1E9-day_smpl)*sampling_rate;
			b_offset_seg=0;
			if(b_offset_data<0)
			{
				b_offset_seg=-b_offset_data;
				b_offset_data=0;
			}

			e_offset_data=((day_smpl+sec_1d)-seg->endtime/1E9)*sampling_rate;
			e_offset_seg=0;
			if(e_offset_data<0)
			{
				e_offset_seg=-e_offset_data;
				e_offset_data=0;
			}
			npts=seg->numsamples-b_offset_seg-e_offset_seg;

			for(i=0;i<npts;i++)
			{
				if(mask[i+b_offset_data]==0)
				{
					//data is empty -> fill with segment
					data[i+b_offset_data]=*((double *)seg->datasamples+i+b_offset_seg);
					mask[i+b_offset_data]=1;
				}
				else if(mask[i+b_offset_data]==1)
				{
					//data already occpuies -> check mismatch between segment and data
					if(data[i+b_offset_data]!=*((double *)seg->datasamples+i+b_offset_seg))
					{
						//mismatch between segment and data
						mask[i+b_offset_data]=2;
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

int read_data(char *stnm,time_t day_smpl,double *data,int *mask)
{
	struct tm day_tm;
	char filepath[50];
	glob_t matches;
	MS3TraceList *mstl=NULL;

	int status;
	int i=0;

	gmtime_r(&day_smpl,&day_tm);

	//construct filepath pattern
	sprintf(filepath,"/home/tllc46/NSNAS/%s/HHZ/NS.%s.HHZ.%d.%03d*",stnm,stnm,1900+day_tm.tm_year,1+day_tm.tm_yday);

	//glob pattern match
	status=glob(filepath,glob_flags,NULL,&matches);
	if(status==GLOB_NOMATCH)
	{
		fprintf(stderr,"No matches at all\n");
		return -1;
	}
	else if(status)
	{
		fprintf(stderr,"Error in glob pattern matching\n");
		return 1;
	}
	/*printf("%ld\n",matches.gl_pathc);
	printf("%s\n",matches.gl_pathv[0]);*/

	//read mseed file to trace list
	while(i<matches.gl_pathc)
	{
		if(ms3_readtracelist(&mstl,matches.gl_pathv[i],NULL,splitversion,mstl_flags,verbose))
		{
			fprintf(stderr,"Error in reading mseed file to trace list\n");
			return 1;
		}
		i++;
	}

	//initialize mask array
	memset(mask,0,npts_1d*sizeof(int));

	//merge segments
	if(0<merge_seg(mstl,day_smpl,data,mask))
	{
		fprintf(stderr,"Error in merging segments\n");
	}

	mstl3_free(&mstl,freeprvtptr);

	return 0;
}

void rtr(double *data,double *detrend)
{
	double denominator=delta*delta*(nfft-1)*nfft*(nfft+1)/12;
	double numerator=0;
	double mean=0;
	double slope;

	int i;

	for(i=0;i<nfft;i++)
	{
		mean+=data[i];
	}
	mean=mean/nfft;

	for(i=0;i<nfft;i++)
	{
		numerator+=deviation[i]*(data[i]-mean);
	}

	slope=numerator/denominator;

	for(i=0;i<nfft;i++)
	{
		detrend[i]=data[i]-mean-slope*deviation[i];
	}
}

void calculate_psd(double *data,fftw_plan *p,double *input,fftw_complex *output,double *psd,FILE *fp)
{
	int i,j;

	memset(psd,0,nfreq_no_zero*sizeof(double));

	for(i=0;i<nwin;i++)
	{
		rtr(data+i*nhop,input);
		fftw_execute(*p);

		for(j=0;j<nfreq_no_zero;j++)
		{
			psd[j]+=output[j+1][0]*output[j+1][0]+output[j+1][1]*output[j+1][1];
		}
	}

	for(i=0;i<nfreq_no_zero;i++)
	{
		psd[i]=psd[i]/nwin;
		psd[i]=log10(psd[i]);
	}

	fwrite(psd,sizeof(double),nfreq_no_zero,fp);
}

void next_day(time_t *cur_day_smpl,struct tm *cur_day_tm,double *data,int *mask)
{
	//shift day
	*cur_day_smpl+=sec_1d;
	gmtime_r(cur_day_smpl,cur_day_tm);

	//shift data and mask array
	memmove(data,data+npts_1d,npts_1d*sizeof(double));
	memmove(mask,mask+npts_1d,npts_1d*sizeof(int));
}

int main(int argc,char **argv)
{
	init_global();
	
	time_t cur_day_smpl;
	struct tm cur_day_tm;
	FILE *fp=fopen("foo1.npy","wb");
	double *data=(double *)malloc(2*npts_1d*sizeof(double));
	int *mask=(int *)malloc(2*npts_1d*sizeof(int));
	double *input=fftw_alloc_real(nfft);
	fftw_complex *output=fftw_alloc_complex(nfreq);
	fftw_plan p=fftw_plan_dft_r2c_1d(nfft,input,output,FFTW_MEASURE);
	double *psd=(double *)malloc(nfreq_no_zero*sizeof(double));

	int status;
	int i,j;

	clock_t x_begin,x_end;
	double cpu_time_used;
	x_begin=clock();

	init_date("2021-11-02","2021-11-02",&cur_day_smpl,&cur_day_tm);

	init_npy(fp);

	//read initial mseed file to trace list
	if(0<read_data(argv[1],cur_day_smpl,data,mask))
	{
		fprintf(stderr,"Error in reading data\n");
		return 1;
	}

	//main loop
	for(i=0;i<nday;i++)
	{
		/*printf("%d %d %d %d %d %d\n",begin_tm.tm_sec,begin_tm.tm_min,begin_tm.tm_hour,begin_tm.tm_mday,begin_tm.tm_mon,begin_tm.tm_year);
		printf("%d %d\n",begin_tm.tm_wday,begin_tm.tm_yday);
		printf("%d %ld %s\n",begin_tm.tm_isdst,begin_tm.tm_gmtoff,begin_tm.tm_zone);
		printf("%ld\n",begin_smpl);*/

		printf("%d-%02d-%02d starting\n",1900+cur_day_tm.tm_year,1+cur_day_tm.tm_mon,cur_day_tm.tm_mday);

		//read next day mseed file to trace list
		if(0<read_data(argv[1],cur_day_smpl+sec_1d,data+npts_1d,mask+npts_1d))
		{
			fprintf(stderr,"Error in reading data\n");
			return 1;
		}

		for(j=0;j<nseg_1d;j++)
		{
			calculate_psd(data+j*nstep,&p,input,output,psd,fp);
		}

		next_day(&cur_day_smpl,&cur_day_tm,data,mask);
	}

	if(fclose(fp))
	{
		fprintf(stderr,"Error in closing npy file\n");
		return 1;
	}
	fftw_destroy_plan(p);
	free(data);
	free(mask);
	free(psd);
	fftw_free(input);
	fftw_free(output);
	term_global();

	x_end=clock();
	cpu_time_used=((double)(x_end-x_begin))/CLOCKS_PER_SEC;
	printf("cpu time used=%f sec\n",cpu_time_used);
	return 0;
}
