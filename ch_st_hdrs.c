/* 컴파일 방법
gcc -o program ch_st_hdrs.c $(bash sac-config -c -l sacio) -lm
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sacio.h> //필수

extern sac *current; //get, set 함수를 쓰지 않고 바로 SAC 구조체 활용 가능

int main(int argc,char **argv)
{
	int max=300000; //npts보다 충분히 큰 수
	int npts,nerr;
	float x[1],y[300000],b,delta,stla,stlo;
	char filer[40]="HL.H35..HHZ.D.2021.276.000000.SAC",filew[8]="foo.sac";
	rsac1(filer,y,&npts,&b,&delta,&max,&nerr,strlen(filer));
	printf("original stla=%.15f\n",current->h->_stla);
	printf("original stlo=%.15f\n",current->h->_stlo);
	printf("original kstnm=%s\n",current->h->kstnm);
	printf("original kevnm=%s\n",current->h->kevnm);	
	if(nerr!=0)
	{
		fprintf(stderr,"Error reading in SAC file: %d\n",nerr);
		return nerr;
	}
	stla=12.345;
	stlo=54.321;
	setfhv("stla",&stla,&nerr,strlen("stla"));
	if(nerr!=0)
	{
		fprintf(stderr,"Error setting stla: %d\n",nerr);
		return nerr;
	}
	setfhv("stlo",&stlo,&nerr,strlen("stlo"));
	if(nerr!=0)
	{
		fprintf(stderr,"Error setting stlo: %d\n",nerr);
		return nerr;
	}

	/* SAC 구조체 직접 활용
	current->h->_stla=12.345;
	current->h->_stlo=54.321;
	*/

	setkhv("kstnm","foo",&nerr,strlen("kstnm"),strlen("foo"));
	if(nerr!=0)
	{
		fprintf(stderr,"Error setting kstnm: %d\n",nerr);
		return nerr;
	}
	setkhv("kevnm","bar",&nerr,strlen("kevnm"),strlen("bar"));
	if(nerr!=0)
	{
		fprintf(stderr,"Error setting kevnm: %d\n",nerr);
		return nerr;
	}
	printf("after stla=%.15f\n",current->h->_stla);
	printf("after stlo=%.15f\n",current->h->_stlo);
	printf("after kstnm=%s\n",current->h->kstnm);

	wsac0(filew,x,y,&nerr,strlen(filew));
	if(nerr!=0)
	{
		fprintf(stderr,"Error writing SAC File: %d\n",nerr);
		return nerr;
	}
	return 0;
}
