#include <stdio.h>
#include <stdlib.h>

int main(int argc,char **argv)
{
	FILE *fp;
	float f[56];
	double f64[22];
	char sf[15];
	int h[22]={0,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,36,35,32,31,54,55};
	int i,n,nvhdr=7;

	for(n=0;n<argc,n++)
	{
		fp=fopen(argv[n+1],"r+b");
		fread(f,sizeof(float),56,fp);
		for(i=0;i<22;i++)
		{
			if(i==0)
			{
				sprintf(sf,"%f",f[h[i]]);
				f64[i]=atof(sf);
			}
			else
			{
				sprintf(sf,"%d",(int)f[h[i]]);
				f64[i]=atof(sf);
			}
		}

		fseek(fp,14*sizeof(float)+6*sizeof(int),SEEK_CUR);
		fwrite(&nvhdr,sizeof(nvhdr),1,fp);
		fseek(fp,0,SEEK_END);
		fwrite(f64,sizeof(double),22,fp);
		fclose(fp);
	}

	return 0;
}
