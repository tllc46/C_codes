#include <stdio.h>
#include <stdlib.h>

int main(int argc,char **argv)
{
	FILE *fp;
	double f[22];
	double delta=0.004;      //delta
	int b=0,e=86400,nvhdr=7; //b,e
	int i,j;
	for(i=1;i<argc;i++)
	{
		fp=fopen(argv[i],"r+b");
		f[0]=delta;
		f[1]=b;
		f[2]=e;
		for(j=3;j<22;j++) f[j]=-12345;
		fseek(fp,70*sizeof(float)+6*sizeof(int),SEEK_SET);
		fwrite(&nvhdr,sizeof(nvhdr),1,fp);
		fseek(fp,0,SEEK_END);
		fwrite(f,sizeof(double),22,fp);
		fclose(fp);
	}
	return 0;
}
