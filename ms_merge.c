/* 컴파일 방법
gcc ms_merge.c -I/usr/local/sac/include -L/usr/local/sac/lib -lm -lmseed -lsacio
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <libmseed.h>
#include <sacio.h>

int main(int argc,char **argv)
{
	MS3TraceList *mstl=NULL;
	MS3TraceID *tid=NULL;
	MS3TraceSeg *seg=NULL;

	int retcode;
	int seg_num=1;
	int b_edge,e_edge;
	int b_offset,e_offset;
	int npts;
	int *num_cover;
	int i;
	int zyear,zjday,zhour,zmin,zsec,zmsec;
	int nerr;

	int8_t splitversion=0;
	int8_t verbose=0;
	int8_t truncate=0;
	int8_t freeprvtptr=1;
	uint8_t hour;
	uint8_t min;
	uint8_t sec;
	uint32_t nsec;
	uint16_t year;
	uint16_t yday;
	uint32_t flags=0;
	int64_t numsamples;

	nstime_t ns_b,ns_e;
        nstime_t ns_b_round,ns_e_round;

	float *data,*xdummy;
	float b=0,delta;
	double sampling_rate=200;

	char msfile1[]="05.SS05.HHZ.2013.304";
	char msfile2[]="05.SS05.HHZ.2013.305";
	char begin[27]="2013-10-31T18:00:00.000000";
	char end[27]="2013-10-31T23:59:59.995000";
	char str_start[27];
	char str_end[27];
	char kname[]="foo.sac";

	ns_b=ms_timestr2nstime(begin);
	ns_e=ms_timestr2nstime(end);

	ns_b_round=((long double)round(ns_b/1E9*sampling_rate)/sampling_rate)*1E9; //소숫점 유효 숫자 유지를 위해 long double로 형변환
        ns_e_round=((long double)round(ns_e/1E9*sampling_rate)/sampling_rate)*1E9;

	if(!(mstl=mstl3_init(NULL)))
	{
		printf("failed to initialize trace list\n");
		return 1;
	}

	flags|=MSF_RECORDLIST;

	retcode=ms3_readtracelist_timewin(&mstl,msfile1,NULL,ns_b_round,ns_e_round,splitversion,flags,verbose);
	if(retcode!=MS_NOTSEED && retcode!=MS_NOERROR)
	{
		printf("failed to read %s\n",msfile1);
		return 1;
	}
	else if(retcode==MS_NOTSEED)
	{
		printf("no time window matches in %s\n",msfile1);
	}

	retcode=ms3_readtracelist_timewin(&mstl,msfile2,NULL,ns_b_round,ns_e_round,splitversion,flags,verbose);
	if(retcode!=MS_NOTSEED && retcode!=MS_NOERROR)
	{
		printf("failed to read %s\n",msfile2);
		return 1;
	}
	else if(retcode==MS_NOTSEED)
	{
		printf("no time window matches in %s\n",msfile2);
	}

	printf("<MS3TraceList>\n");
	printf("numtraceids=%d\n",mstl->numtraceids);
	printf("\n");

	tid=mstl->traces.next[0];
	while(tid)
	{
		tid=mstl->traces.next[0];
		if(!ms_nstime2timestr(tid->earliest,str_start,ISOMONTHDAY,MICRO))
		{
			printf("failed to convert trace ID earliest\n");
			return 1;
		}
		if(!ms_nstime2timestr(tid->latest,str_end,ISOMONTHDAY,MICRO))
		{
			printf("failed to convert trace ID latest\n");
			return 1;
		}

		printf("<MS3TraceID>\n");
		printf("sid=%s\n",tid->sid);
		printf("pubversion=%d\n",tid->pubversion);
		printf("earliest=%s\n",str_start);
		printf("latest=%s\n",str_end);
		printf("numsegments=%d\n",tid->numsegments);
		printf("height=%d\n",tid->height);
		printf("\n");

		seg=tid->first;
		printf("<MS3TraceSeg>\n");
		while(seg)
		{
			printf("<MS3TraceSeg %d>\n",seg_num);
			if(!ms_nstime2timestr(seg->starttime,str_start,ISOMONTHDAY,MICRO))
			{
				printf("failed to convert trace segment starttime\n");
				return 1;
			}
			if(!ms_nstime2timestr(seg->endtime,str_end,ISOMONTHDAY,MICRO))
			{
				printf("failed to convert latest to time string\n");
				return 1;
			}

			printf("starttime=%s\n",str_start);
			printf("endtime=%s\n",str_end);
			printf("samprate=%f\n",seg->samprate);
			printf("samplecnt=%ld\n",seg->samplecnt);
			printf("\n");

			seg=seg->next;
			seg_num++;
		}
		tid=tid->next[0];
	}

	npts=(ns_e_round-ns_b_round)/1E9*sampling_rate+1; //마지막 시간 포함하므로 +1

	data=(float *)malloc(sizeof(float)*npts);
	num_cover=(int *)malloc(sizeof(int)*npts);

	for(i=0;i<npts;i++)
	{
		num_cover[i]=0;
	}

	tid=mstl->traces.next[0];
	while(tid)
	{
		seg=tid->first;
		while(seg)
		{
			if(mstl3_unpack_recordlist(tid,seg,NULL,0,verbose)==-1)
			{
				printf("fail to unpack record list\n");
				return 1;
			}

			if(mstl3_convertsamples(seg,'f',truncate)!=MS_NOERROR)
			{
				printf("fail converting data samples from integer to float\n");
				return 1;
			}

			b_edge=(seg->starttime-ns_b_round)*sampling_rate/1E9;
			e_edge=(ns_e_round-seg->endtime)*sampling_rate/1E9;
			b_offset=0;
			e_offset=0;
			if(b_edge<0)
			{
				b_offset=-b_edge;
				b_edge=0;
			}
			if(e_edge<0)
			{
				e_offset=-e_edge;
				e_edge=0;
			}

			numsamples=seg->numsamples-b_offset-e_offset;
			memcpy(data+b_edge,(float *)seg->datasamples+b_offset,sizeof(float)*numsamples);
			for(i=b_edge;i<b_edge+numsamples;i++)
			{
				num_cover[i]++;
			}

			seg=seg->next;
		}
		tid=tid->next[0];
	}

	for(i=0;i<npts;i++)
	{
		if(!num_cover[i])
		{
			data[i]=-99;
		}
	}

	delta=1/sampling_rate;
	if(ms_nstime2time(ns_b_round,&year,&yday,&hour,&min,&sec,&nsec)!=MS_NOERROR)
	{
		printf("fail to convert nanoseconds time to time\n");
		return 1;
	}
	zyear=year;
	zjday=yday;
	zhour=hour;
	zmin=min;
	zsec=sec;
	zmsec=nsec/1E6;

	newhdr();
	setnhv("npts",&npts,&nerr,strlen("npts"));
	setnhv("nzyear",&zyear,&nerr,strlen("nzyear"));
	setnhv("nzjday",&zjday,&nerr,strlen("nzjday"));
	setnhv("nzhour",&zhour,&nerr,strlen("nzhour"));
	setnhv("nzmin",&zmin,&nerr,strlen("nzmin"));
	setnhv("nzsec",&zsec,&nerr,strlen("nzsec"));
	setnhv("nzmsec",&zmsec,&nerr,strlen("nzmsec"));
	setfhv("b",&b,&nerr,strlen("b"));
	setfhv("delta",&delta,&nerr,strlen("delta"));
	wsac0(kname,xdummy,data,&nerr,strlen(kname));

	if(nerr)
	{
		printf("fail writing SAC file\n");
		return 1;
	}

	mstl3_free(&mstl,freeprvtptr);
	free(data);
	free(num_cover);

	return 0;
}
