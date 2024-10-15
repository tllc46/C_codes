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
int rfft_scale=2;
double dtiny=1e-21;

int sampling_rate,avg_len;
double avg_ovrlp,freq_smoothing_width_octaves,freq_step_octaves;
char begin_str[11],end_str[11],data_path[128],resp_path[128],output_path[128];

int npts_1d;
double delta;

double *data;
int *mask;

int avg_npts,avg_shift_npts,navg_1d;

int nfft,sub_shift_npts,nsub,nfreq,nfreq_no_zero;

double *fftw_input;
fftw_complex *fftw_output;
fftw_plan plan;

double *deviation;

int ntaper;
double *cos_taper,taper_scale;

int nfreq_bin;
double *spectrum,*psd,*psd_gap;

int *right_idx,*left_idx;

double *resp_amp;

time_t cur_day_sct;
struct tm cur_day_bdt;
int nday,navg;

FILE *fp_out;

int cur_day_gap,next_day_gap;

int glob_flags=0;
int8_t splitversion=0;
uint32_t mstl_flags=0;
int8_t verbose=0;
int8_t truncate=0;
int8_t freeprvtptr=1;

void read_global(void)
{
	char line[128];
	FILE *fp_in;

	fp_in=fopen("ppsd.in","r");

	fgets(line,128,fp_in);
	sscanf(line,"%d",&sampling_rate);

	fgets(line,128,fp_in);
	sscanf(line,"%d",&avg_len);

	fgets(line,128,fp_in);
	sscanf(line,"%lf",&avg_ovrlp);

	fgets(line,128,fp_in);
	sscanf(line,"%lf",&freq_smoothing_width_octaves);

	fgets(line,128,fp_in);
	sscanf(line,"%lf",&freq_step_octaves);

	fgets(line,128,fp_in);
	sscanf(line,"%s",begin_str);

	fgets(line,128,fp_in);
	sscanf(line,"%s",end_str);

	fgets(line,128,fp_in);
	sscanf(line,"%s",data_path);

	fgets(line,128,fp_in);
	sscanf(line,"%s",resp_path);

	fgets(line,128,fp_in);
	sscanf(line,"%s",output_path);

	fclose(fp_in);
}

void init_global(void)
{
	int avg_shift,win_len,exponent;

	int i;

	npts_1d=sec_1d*sampling_rate;
	delta=(double)1/sampling_rate;

	//initialize data and mask arrays
	data=(double *)malloc(2*npts_1d*sizeof(double));
	mask=(int *)malloc(2*npts_1d*sizeof(int));

	//initialize average window
	avg_npts=avg_len*sampling_rate;
	avg_shift=avg_len*(1-avg_ovrlp);
	avg_shift_npts=avg_shift*sampling_rate;
	navg_1d=sec_1d/avg_shift;

	//initialize sub window
	//13 sub windows overlapping by 75%
	//1 sub window + 25% * 12 sub windows
	win_len=avg_len*sampling_rate/4;
	exponent=log2(win_len);
	nfft=exp2(exponent); //previous power of 2
	sub_shift_npts=nfft/4; //75% overlap -> 25% shift
	nsub=(avg_npts-nfft)/sub_shift_npts+1;
	
	//initialize FFT
	nfreq=nfft/2+1;
	nfreq_no_zero=nfreq-1;
	fftw_input=fftw_alloc_real(nfft);
	fftw_output=fftw_alloc_complex(nfreq);
	plan=fftw_plan_dft_r2c_1d(nfft,fftw_input,fftw_output,FFTW_MEASURE);

	//initialize detrend
	deviation=(double *)malloc(nfft*sizeof(double));
	for(i=0;i<nfft/2;i++)
	{
		deviation[nfft/2+i]=0.5*(2*i+1)*delta;
		deviation[(nfft/2-1)-i]=-deviation[nfft/2+i];
	}

	//initialize cosine taper
	ntaper=0.1*nfft;
	cos_taper=(double *)malloc(ntaper*sizeof(double));
	for(i=0;i<ntaper;i++)
	{
		cos_taper[i]=0.5*(1-cos(M_PI*i/(ntaper-1)));
	}
	taper_scale=0.875*nfft;

	//initialize psd results
	nfreq_bin=(exponent-1)/freq_step_octaves+1;
	spectrum=(double *)malloc(nfreq_no_zero*sizeof(double));
	psd=(double *)malloc(nfreq_bin*sizeof(double));
	psd_gap=(double *)malloc(nfreq_bin*sizeof(double));
	for(i=0;i<nfreq_bin;i++)
	{
		psd_gap[i]=-999999;
	}

	//initialize read trace
	mstl_flags|=MSF_UNPACKDATA;
}

void term_global(void)
{
	free(data);
	free(mask);

	fftw_free(fftw_input);
	fftw_free(fftw_output);
	fftw_destroy_plan(plan);

	free(deviation);
	
	free(cos_taper);

	free(spectrum);
	free(psd);
	free(psd_gap);

	free(right_idx);
	free(left_idx);

	free(resp_amp);

	fclose(fp_out);
}

void init_freq_smooth(void)
{
	double *log2_freq;
	double right_exp,left_exp;

	int i,jr=0,jl=0;

	//log2(frequency)
	log2_freq=(double *)malloc(nfreq_no_zero*sizeof(double));
	for(i=0;i<nfreq_no_zero;i++)
	{
		log2_freq[i]=log2(i+1);
	}

	//indices of spectral smoothing range in frequency
	right_idx=(int *)malloc(nfreq_bin*sizeof(int));
	left_idx=(int *)malloc(nfreq_bin*sizeof(int));
	for(i=0;i<nfreq_bin;i++)
	{
		right_exp=i*freq_step_octaves+0.5*freq_smoothing_width_octaves;
		while(log2_freq[jr]<=right_exp && jr<nfreq_no_zero)
		{
			jr++;
		}
		right_idx[i]=jr;

		left_exp=i*freq_step_octaves-0.5*freq_smoothing_width_octaves;
		while(log2_freq[jl]<left_exp && jl<nfreq_no_zero)
		{
			jl++;
		}
		left_idx[i]=jl;
	}

	free(log2_freq);
}

void init_response(void)
{
	evalresp_options *options=NULL;
	evalresp_filter *filter=NULL;
	evalresp_channels *channels=NULL;
	evalresp_response *response=NULL;

	int i;

	//initialize evalresp options
	evalresp_new_options(NULL,&options);
	evalresp_new_filter(NULL,&filter);
	options->min_freq=0;
	options->max_freq=sampling_rate/2;
	options->nfreq=nfreq;
	options->lin_freq=1;
	options->format=evalresp_complex_output_format;
	options->unit=evalresp_acceleration_unit;

	//calculate instrumental response
	evalresp_filename_to_channels(NULL,resp_path,options,filter,&channels);
	evalresp_channel_to_response(NULL,channels->channels[0],options,&response);

	//calculate amplitude of instrumental response
	resp_amp=(double *)malloc(nfreq_no_zero*sizeof(double));
	for(i=0;i<nfreq_no_zero;i++)
	{
		resp_amp[i]=response->rvec[1+i].real*response->rvec[1+i].real+response->rvec[1+i].imag*response->rvec[1+i].imag;
	}

	evalresp_free_options(&options);
	evalresp_free_filter(&filter);
	evalresp_free_channels(&channels);
	evalresp_free_response(&response);
}

void init_date()
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
}

void init_npy(void)
{
	char magic[7]="\x93" "NUMPY";
	char major_ver=1,minor_ver=0;
	unsigned short header_len=118;
	char *header=(char *)malloc((header_len+1)*sizeof(char));
	int status;

	fp_out=fopen(output_path,"wb");

	memset(header,' ',header_len);

	status=sprintf(header,"{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d), }",navg,nfreq_bin);

	header[status]=' ';
	header[header_len-1]='\n';
	header[header_len]='\0';

	fputs(magic,fp_out);

	fwrite(&major_ver,sizeof(char),1,fp_out);
	fwrite(&minor_ver,sizeof(char),1,fp_out);
	fwrite(&header_len,sizeof(unsigned short),1,fp_out);

	fputs(header,fp_out);
}

int merge_seg(MS3TraceList *mstl,time_t cur_day_sct,double *data,int *mask)
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

			data_loff=round((seg->starttime/1E9-cur_day_sct)*sampling_rate);
			seg_loff=0;
			if(data_loff<0)
			{
				//segment left part lies outside data
				seg_loff=-data_loff;
				data_loff=0;
			}

			data_roff=round(((cur_day_sct+sec_1d)-seg->endtime/1E9)*sampling_rate)-1;
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

int read_mseed(time_t cur_day_sct,double *data,int *mask)
{
	struct tm cur_day_bdt;
	char pattern[68];
	glob_t matches;
	MS3TraceList *mstl=NULL;

	int status;
	int i=0;

	gmtime_r(&cur_day_sct,&cur_day_bdt);

	//construct pattern
	sprintf(pattern,"%s.%d.%03d*",data_path,1900+cur_day_bdt.tm_year,1+cur_day_bdt.tm_yday);

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
	memset(mask,0,npts_1d*sizeof(int));

	//merge segments
	if(0<merge_seg(mstl,cur_day_sct,data,mask))
	{
		fprintf(stderr,"Error in merging segments\n");
		return -1;
	}

	mstl3_free(&mstl,freeprvtptr);

	return 0;
}

void detrend(double *data)
{
	double denominator=delta*delta*(nfft-1)*nfft*(nfft+1)/12;
	double numerator=0;
	double mean=0;
	double slope;

	int i;

	//mean_y
	for(i=0;i<nfft;i++)
	{
		mean+=data[i];
	}
	mean=mean/nfft;

	//covariance_xy
	for(i=0;i<nfft;i++)
	{
		numerator+=deviation[i]*(data[i]-mean);
	}

	slope=numerator/denominator;

	for(i=0;i<nfft;i++)
	{
		fftw_input[i]=data[i]-mean-slope*deviation[i];
	}
}

void taper(double *data)
{
	int i;

	for(i=0;i<ntaper;i++)
	{
		fftw_input[i]=data[i]*cos_taper[i];
		fftw_input[(nfft-1)-i]=data[(nfft-1)-i]*cos_taper[i];
	}
}

void calculate_psd(void)
{
	int gap_flag;
	int i,j,k;

	for(i=0;i<navg_1d;i++)
	{
		gap_flag=0;
		for(j=0;j<avg_npts;j++)
		{
			if(mask[i*avg_shift_npts+j]!=1)
			{
				//average window contains gap
				gap_flag=1;
				break;
			}
		}

		if(gap_flag)
		{
			fprintf(stderr,"%02d/%d average window gap exists\n",i+1,navg_1d);
			fwrite(psd_gap,sizeof(double),nfreq_bin,fp_out);
			continue;
		}

		memset(spectrum,0,nfreq_no_zero*sizeof(double));
		memset(psd,0,nfreq_bin*sizeof(double));

		for(j=0;j<nsub;j++)
		{
			detrend(data+i*avg_shift_npts+j*sub_shift_npts);
			taper(fftw_input);
			fftw_execute(plan);

			for(k=0;k<nfreq_no_zero;k++)
			{
				spectrum[k]+=fftw_output[1+k][0]*fftw_output[1+k][0]+fftw_output[1+k][1]*fftw_output[1+k][1];
			}
		}

		for(j=0;j<nfreq_no_zero;j++)
		{
			spectrum[j]=spectrum[j]*delta/(nsub*taper_scale*resp_amp[j]);
			if(j!=nfreq_no_zero-1)
			{
				spectrum[j]=rfft_scale*spectrum[j];
			}
			if(spectrum[j]<dtiny)
			{
				spectrum[j]=dtiny;
			}
			spectrum[j]=10*log10(spectrum[j]);
		}

		for(j=0;j<nfreq_bin;j++)
		{
			for(k=left_idx[j];k<right_idx[j];k++)
			{
				psd[j]+=spectrum[k];
			}
			psd[j]=psd[j]/(right_idx[j]-left_idx[j]);
		}

		fwrite(psd,sizeof(double),nfreq_bin,fp_out);
	}
}

void next_day(void)
{
	//shift day
	cur_day_sct+=sec_1d;
	gmtime_r(&cur_day_sct,&cur_day_bdt);

	cur_day_gap=next_day_gap;

	if(!next_day_gap)
	{
		//shift data and mask array
		memmove(data,data+npts_1d,npts_1d*sizeof(double));
		memmove(mask,mask+npts_1d,npts_1d*sizeof(int));
	}
}

int main(int argc,char **argv)
{
	char info[20];

	int i,j;

	init_global();
	init_freq_smooth();
	init_response();
	init_date();
	init_npy();

	//read initial mseed file to trace list
	if((cur_day_gap=read_mseed(cur_day_sct,data,mask))==-1)
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
		if((next_day_gap=read_mseed(cur_day_sct+sec_1d,data+npts_1d,mask+npts_1d))==-1)
		{
			fprintf(stderr,"Error in reading data\n");
			return 1;
		}

		if(cur_day_gap==1)
		{
			fprintf(stderr,"no file at all\n");
			for(j=0;j<navg_1d;j++)
			{
				fwrite(psd_gap,sizeof(double),nfreq_bin,fp_out);
			}
			next_day();
			continue;
		}

		calculate_psd();

		next_day();
	}

	term_global();

	return 0;
}
