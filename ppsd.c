/*
gcc ppsd.c -lmseed -lfftw3 -levalresp -levalresp_log -lspline -lmxmlev -lm
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
#include <evalresp/public_api.h>

int sec_1d=86400;
int sampling_rate=200;
int ppsd_length=3600;
double overlap=0.5;
int period_smoothing_width_octaves=1;
double period_step_octaves=0.125;

int npts_1d;
double delta;
int seg_len,nseg_1d,nstep;
int nfft,nhop,nwin,nfreq,ntaper;
int nfreq_no_zero,nfreq_bin;
int *right_idx,*left_idx;
double *deviation;
double *cos_taper,taper_factor=0;
double *respamp;
double *data;
int *mask;
double *input;
fftw_complex *output;
fftw_plan p;
double *spec,*psd,*gap_psd;
int scaling_factor=2;

time_t cur_day_smpl;
struct tm cur_day_tm;
int nday,nseg;

FILE *fp;

int glob_flags=0;
int8_t splitversion=0;
uint32_t mstl_flags=0;
int8_t verbose=0;
int8_t truncate=0;
int8_t freeprvtptr=1;

void init_global(void)
{
	int step,win_len,exponent;

	int i;

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
	nfreq_no_zero=nfreq-1;
	nfreq_bin=(exponent-1)/period_step_octaves+1;
	ntaper=0.1*nfft;

	//initialize deviation
	deviation=(double *)malloc(nfft*sizeof(double));
	for(i=0;i<nfft/2;i++)
	{
		deviation[i+nfft/2]=0.5*(2*i+1)*delta;
		deviation[-i+(nfft/2-1)]=-deviation[i+nfft/2];
	}

	//initialize cosine taper function
	cos_taper=(double *)malloc(ntaper*sizeof(double));
	for(i=0;i<ntaper;i++)
	{
		cos_taper[i]=0.5*(1-cos(M_PI*i/(ntaper-1)));
		taper_factor+=cos_taper[i]*cos_taper[i];
	}
	taper_factor=2*taper_factor;
	taper_factor+=nfft-2*ntaper;

	//data and mask arrays
	data=(double *)malloc(2*npts_1d*sizeof(double));
	mask=(int *)malloc(2*npts_1d*sizeof(int));

	//ffts
	input=fftw_alloc_real(nfft);
	output=fftw_alloc_complex(nfreq);
	p=fftw_plan_dft_r2c_1d(nfft,input,output,FFTW_MEASURE);

	//psd results
	spec=(double *)malloc(nfreq_no_zero*sizeof(double));
	psd=(double *)malloc(nfreq_bin*sizeof(double));
	gap_psd=(double *)malloc(nfreq_bin*sizeof(double));
	for(i=0;i<nfreq_bin;i++)
	{
		gap_psd[i]=999999;
	}

	mstl_flags|=MSF_UNPACKDATA;
}

void term_global(void)
{
	free(right_idx);
	free(left_idx);
	free(deviation);
	free(data);
	free(mask);
	fclose(fp);
	fftw_free(input);
	fftw_free(output);
	fftw_destroy_plan(p);
	free(spec);
	free(psd);
	free(gap_psd);
}

void init_smooth(void)
{
	double *frequency_log2;
	double right_exp,left_exp;

	int i,jr=0,jl=0;

	frequency_log2=(double *)malloc(nfreq_no_zero*sizeof(double));
	for(i=0;i<nfreq_no_zero;i++)
	{
		frequency_log2[i]=log2(i+1);
	}

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

	free(frequency_log2);
}

void init_resp(void)
{
	evalresp_options *options=NULL;
	evalresp_filter *filter=NULL;
	evalresp_channels *channels=NULL;
	evalresp_response *response=NULL;

	int i;

	evalresp_new_options(NULL,&options);
	evalresp_new_filter(NULL,&filter);
	options->min_freq=0;
	options->max_freq=sampling_rate/2;
	options->nfreq=nfreq;
	options->lin_freq=1;
	options->format=evalresp_complex_output_format;
	options->unit=evalresp_acceleration_unit;

	evalresp_filename_to_channels(NULL,"/home/tllc46/48NAS1/tllc46/Jeju/RESP/Broadband/RESP.05.SS10..HHZ",options,filter,&channels);
	evalresp_channel_to_response(NULL,channels->channels[0],options,&response);

	respamp=(double *)malloc(nfreq_no_zero*sizeof(double));
	for(i=0;i<nfreq_no_zero;i++)
	{
		respamp[i]=response->rvec[i+1].real*response->rvec[i+1].real+response->rvec[i+1].imag*response->rvec[i+1].imag;
	}

	evalresp_free_options(&options);
	evalresp_free_filter(&filter);
	evalresp_free_channels(&channels);
	evalresp_free_response(&response);
}

void init_date(char *begin_str,char *end_str)
{
	struct tm begin_tm,end_tm;
	time_t begin_smpl,end_smpl;

	strptime(begin_str,"%Y-%m-%d",&begin_tm);
	begin_tm.tm_hour=0;
	begin_tm.tm_min=0;
	begin_tm.tm_sec=0;
	cur_day_tm=begin_tm;

	strptime(end_str,"%Y-%m-%d",&end_tm);
	end_tm.tm_hour=0;
	end_tm.tm_min=0;
	end_tm.tm_sec=0;

	begin_smpl=timegm(&begin_tm);
	cur_day_smpl=begin_smpl;
	end_smpl=timegm(&end_tm);
	nday=(end_smpl-begin_smpl)/sec_1d+1;
	nseg=nday*nseg_1d;
}

void init_npy(void)
{
	char magic[7]="\x93" "NUMPY";
	char major=1,minor=0;
	unsigned short header_len=118;
	char *header=(char *)malloc((header_len+1)*sizeof(char));
	int status;

	fp=fopen("foo1.npy","wb");

	memset(header,' ',header_len);

	status=sprintf(header,"{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d), }",nseg,nfreq_bin);

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

			b_offset_data=round((seg->starttime/1E9-day_smpl)*sampling_rate);
			b_offset_seg=0;
			if(b_offset_data<0)
			{
				b_offset_seg=-b_offset_data;
				b_offset_data=0;
			}

			e_offset_data=round(((day_smpl+sec_1d)-seg->endtime/1E9)*sampling_rate);
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
	char filepath[68];
	glob_t matches;
	MS3TraceList *mstl=NULL;

	int status;
	int i=0;

	gmtime_r(&day_smpl,&day_tm);

	//construct filepath pattern
	sprintf(filepath,"/home/tllc46/48NAS2/symbolic.jeju/05.%s/HHZ/05.%s.HHZ.%d.%03d*",stnm,stnm,1900+day_tm.tm_year,1+day_tm.tm_yday);

	//glob pattern match
	status=glob(filepath,glob_flags,NULL,&matches);
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
	memset(mask,0,npts_1d*sizeof(int));

	//merge segments
	if(0<merge_seg(mstl,day_smpl,data,mask))
	{
		fprintf(stderr,"Error in merging segments\n");
		return -1;
	}

	mstl3_free(&mstl,freeprvtptr);

	return 0;
}

void rtr(double *data)
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
		input[i]=data[i]-mean-slope*deviation[i];
	}
}

void taper(double *data)
{
	int i;

	for(i=0;i<ntaper;i++)
	{
		input[i]=data[i]*cos_taper[i];
		input[(nfft-1)-i]=data[(nfft-1)-i]*cos_taper[i];
	}
}

void calculate_psd(void)
{
	int gap_flag;
	int i,j,k;

	for(i=0;i<nseg_1d;i++)
	{
		gap_flag=0;
		for(j=0;j<seg_len;j++)
		{
			if(mask[i*nstep+j]!=1)
			{
				gap_flag=1;
				break;
			}
		}

		if(gap_flag)
		{
			fprintf(stderr,"%02d/%d segment gap exists\n",i+1,nseg_1d);
			fwrite(gap_psd,sizeof(double),nfreq_bin,fp);
			continue;
		}

		memset(spec,0,nfreq_no_zero*sizeof(double));
		memset(psd,0,nfreq_bin*sizeof(double));

		for(j=0;j<nwin;j++)
		{
			rtr(data+i*nstep+j*nhop);
			taper(input);
			fftw_execute(p);

			for(k=0;k<nfreq_no_zero;k++)
			{
				spec[k]+=output[k+1][0]*output[k+1][0]+output[k+1][1]*output[k+1][1];
			}
		}

		for(j=0;j<nfreq_no_zero;j++)
		{
			spec[j]=spec[j]*delta/(nwin*taper_factor*respamp[j]);
			if(j!=nfreq_no_zero-1)
			{
				spec[j]=scaling_factor*spec[j];
			}
			spec[j]=10*log10(spec[j]);
		}

		for(j=0;j<nfreq_bin;j++)
		{
			for(k=left_idx[j];k<right_idx[j];k++)
			{
				psd[j]+=spec[k];
			}
			psd[j]=psd[j]/(right_idx[j]-left_idx[j]);
		}

		fwrite(psd,sizeof(double),nfreq_bin,fp);
	}
}

void next_day(int *cur_day_gap,int *next_day_gap)
{
	//shift day
	cur_day_smpl+=sec_1d;
	gmtime_r(&cur_day_smpl,&cur_day_tm);

	*cur_day_gap=*next_day_gap;

	if(!*next_day_gap)
	{
		//shift data and mask array
		memmove(data,data+npts_1d,npts_1d*sizeof(double));
		memmove(mask,mask+npts_1d,npts_1d*sizeof(int));
	}
}

int main(int argc,char **argv)
{
	int cur_day_gap,next_day_gap;
	char info[20];

	int i,j;

	init_global();
	init_smooth();
	init_resp();
	init_date("2013-10-01","2015-10-31");
	init_npy();

	//read initial mseed file to trace list
	if((cur_day_gap=read_data(argv[1],cur_day_smpl,data,mask))==-1)
	{
		fprintf(stderr,"Error in reading data\n");
		return 1;
	}

	//main loop
	for(i=0;i<nday;i++)
	{
		strftime(info,20,"%Y-%m-%d starting",&cur_day_tm);
		printf("%s\n",info);

		//read next day mseed file to trace list
		if((next_day_gap=read_data(argv[1],cur_day_smpl+sec_1d,data+npts_1d,mask+npts_1d))==-1)
		{
			fprintf(stderr,"Error in reading data\n");
			return 1;
		}

		if(cur_day_gap==1)
		{
			fprintf(stderr,"no file at all\n");
			for(j=0;j<nseg_1d;j++)
			{
				fwrite(gap_psd,sizeof(double),nfreq_bin,fp);
			}
			next_day(&cur_day_gap,&next_day_gap);
			continue;
		}

		calculate_psd();

		next_day(&cur_day_gap,&next_day_gap);
	}

	term_global();

	return 0;
}
