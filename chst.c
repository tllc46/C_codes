/* 
sac/aux/external/flipxy.c 참고

<방법>
(1)컴파일
gcc -fPIC -I/usr/local/sac/aux/external -c chst.c
gcc -shared  -o libchst.so chst.o /usr/local/sac/lib/libsac_all.a

(2)
~/.bashrc 마지막 줄에
export SACSOLIST=libchst.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/libchst.so
추가

(3) sac에서 사용
sac
load chst
chst arg1 arg2 arg3 ...

(4)
"Storing data without a filename"이 뜨는데 이는 dfm/updatedfl.c/updatedfl()에서
호출한 sacput(s)가 원인
bs/iniam.c/sacput() 원형을 고치면 된다
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extfunc.h" //필수

int chst(int argc,char **argv,sac_files *call_data,int *update)
{
	int i,numfiles,error;
	float stla,stlo;
	char *kstnm;
	sac_header *hdr;
	*update=REPLACE;
	numfiles=call_data->nfiles;
	for(i=0;i<numfiles;i++)
	{
		hdr=call_data->ext_hdrs[i];
		stla=getfhdr(hdr,"stla",&error);
		if(error!=0)
		{
			fprintf(stderr,"Error getting stla: %d\n",error);
			return error;
		}
		printf("initial stla is %f degrees\n",stla);
		stlo=getfhdr(hdr,"stlo",&error);
		if(error!=0)
		{
			fprintf(stderr,"Error getting stlo: %d\n",error);
			return error;
		}
		printf("initial stlo is %f degrees\n",stlo);
		kstnm=getahdr(hdr,"kstnm",&error);
		if(error!=0)
		{
			fprintf(stderr,"Error getting kstnm: %d\n",error);
			return error;
		}
		printf("initial kstnm is %s\n",kstnm);

		stla=atof(argv[1]);
		stlo=atof(argv[2]);

		setfhdr(hdr,"stla",stla,&error);
		if(error!=0)
		{
			fprintf(stderr,"Error setting stla: %d\n",error);
			return error;
		}
		setfhdr(hdr,"stlo",stlo,&error);
		if(error!=0)
		{
			fprintf(stderr,"Error setting stlo: %d\n",error);
			return error;
		}
		setahdr(hdr,"kstnm",argv[3],&error);
		if(error!=0)
		{
			fprintf(stderr,"Error getting kstnm: %d\n",error);
			return error;
		}
	}
	return 0;
}
