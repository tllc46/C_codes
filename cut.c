/*
how to compile
gcc -o program cut.c $(bash sac-config -c -l sac sacio) -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sacio.h> //must include
#include <sac.h> //must include

extern sac *current; //we can use sac struct directly, skipping 
int main(int argc,char **argv)
{
	int max=400000;
	int npts,nerr,nout=max;
	float x[1],y[max],out[max];
	float b,delta;
	char filer[40]="HL.H35..HHZ.D.2021.276.000000.SAC",filew[10]="foo.sac";
	rsac1(filer,y,&npts,&b,&delta,&max,&nerr,strlen(filer));

	cut(y,npts,b,delta,21500,21600,1,out,&nout);

	printf("npts=%d\n",current->h->npts);
	printf("b=%f\n",current->h->_b);
	printf("e=%f\n",current->h->_e);
	printf("nout=%d\n",nout);
	printf("length=%ld\n",sizeof(out));

	current->h->npts=nout;
	current->z->_b=21500.0;
	wsac0(filew,x,out,&nerr,strlen(filew));

	return 0;
}
