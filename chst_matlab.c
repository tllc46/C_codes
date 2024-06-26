//unfinished code!!

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "engine.h" //must include
#include "matrix.h" //must include

int main(int argc,char **argv)
{
	int i,numfiles,error;
	float stla=12.345;
	float stlo=54.321;
	Engine *ep; //matlab engine structure
	mxArray *STLA,*STLO; //matlab array structure
	/*
	mwSize: same as size_t
	mxClassID class id: mxINT32_CLASS, mxSINGLE_CLASS, mxDOUBLE_CLASS, ... etc. 
	mxComplexity ComplexFlag: mxREAL or mxCOMPLEX
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
