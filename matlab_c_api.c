/*
컴파일 방법
(1)
mex -R2018a -client engine matlab_c_api.c

(2)
~/.bashrc 파일 마지막 줄에
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:/usr/local/MATLAB/R2023a/bin/glnxa64:/usr/local/MATLAB/R2023a/sys/os/glnxa64
추가
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "engine.h" //필수
#include "matrix.h" //필수

int main(int argc,char **argv)
{
	int i,numfiles,error;
	float stla=12.345;
	float stlo=54.321;
	Engine *ep; //MatLab engine 구조체
	mxArray *STLA,*STLO; //MatLab array 구조체
	/*
 	mxCreateNumericMatrix(mwSize,mwSize,mxClassID,mxComplexity)
	mwSize: size_t 크기 만큼의 요소 개수
 	mxClassID: mxINT32_CLASS, mxSINGLE_CLASS, mxDOUBLE_CLASS, ... 등
	mxComplexity ComplexFlag: mxREAL 또는 mxCOMPLEX
	*/
	STLA=mxCreateNumericMatrix(1,1,mxSINGLE_CLASS,mxREAL);
	STLO=mxCreateNumericMatrix(1,1,mxSINGLE_CLASS,mxREAL);

	if(!(ep=engOpen("matlab -nodisplay")))
	{
		fprintf(stderr,"Can't open MATLAB engine\n");
	}

	printf("initial stla is %f degrees\n",stla);
	printf("initial stlo is %f degrees\n",stlo);

	memcpy(mxGetSingles(STLA),&stla,sizeof(stla));
	memcpy(mxGetSingles(STLO),&stlo,sizeof(stlo));

	error=engPutVariable(ep,"STLA",STLA);
	if(error)
	{
		fprintf(stderr,"Can't put variable STLA\n");
	}

	error=engPutVariable(ep,"STLO",STLO);
	if(error)
	{
		fprintf(stderr,"Can't put variable STLO\n");
	}
	error=engEvalString(ep,"STLA=2*STLA");
	if(error)
	{
		fprintf(stderr,"Can't evaluate stla expression\n");
	}
	error=engEvalString(ep,"STLO=2*STLO");
	if(error)
	{
		fprintf(stderr,"Can't evaluate stlo expression\n");
	}

	STLA=engGetVariable(ep,"STLA");
	STLO=engGetVariable(ep,"STLO");

	memcpy(&stla,mxGetSingles(STLA),sizeof(stla));
	memcpy(&stlo,mxGetSingles(STLO),sizeof(stlo));

	printf("final stla=%f\n",stlo);
	printf("final stlo=%f\n",stla);

	mxDestroyArray(STLA);
	mxDestroyArray(STLO);
	error=engClose(ep);
	if(error)
	{
		fprintf(stderr,"Can't close MATLAB engine\n");
	}

	return 0;
}
