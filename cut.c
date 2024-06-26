/*
컴파일 방법
gcc -o program cut.c $(bash sac-config -c -l sac sacio) -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sacio.h> //필수
#include <sac.h> //필수

extern sac *current; //get, set 함수를 쓰지 않고 바로 SAC 구조체 활용 가능
int main(int argc,char **argv)
{
	int max=400000; //npts보다 충분히 큰 수
	int npts,nerr,nout=max; //nout 또한 npts보다 충분히 큰 수
	float x[1],y[max],out[max];
	float b,delta;
	char filer[40]="HL.H35..HHZ.D.2021.276.000000.SAC",filew[10]="foo.sac";
	rsac1(filer,y,&npts,&b,&delta,&max,&nerr,strlen(filer)); //sac IO 함수

	cut(y,npts,b,delta,21500,21600,1,out,&nout); //SAC library 함수

	printf("npts=%d\n",current->h->npts);
	printf("b=%f\n",current->h->_b);
	printf("e=%f\n",current->h->_e);
	printf("nout=%d\n",nout);
	printf("length=%ld\n",sizeof(out));

	current->h->npts=nout; //sac 구조체 활용
	current->z->_b=21500.0; //sac 구조체 활용
	wsac0(filew,x,out,&nerr,strlen(filew)); //SAC IO 함수

	return 0;
}
