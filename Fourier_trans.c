/* 컴파일 방법
gcc test3.c -Ifftw/include -I/usr/local/sac/include -Lfftw/lib -L/usr/local/sac/lib -lfftw3f -lsacio -lm
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sacio.h>
#include <fftw3.h>

int main(int argc,char **argv)
{
        int max=300000;
        int npts,nerr,i;
        float x[1],y[300000],b,delta,df;
        char filer[11]="seismo.sac",filew[8]="foo.sac";

        fftwf_plan p;

        rsac1(filer,y,&npts,&b,&delta,&max,&nerr,strlen(filer));

        df=1/(delta*npts);

        float *in=fftwf_alloc_real(npts);
        fftwf_complex *out=fftwf_alloc_complex(npts/2+1);

        p=fftwf_plan_dft_r2c_1d(npts,in,out,FFTW_ESTIMATE);

        for(i=0;i<npts;i++)
        {
                in[i]=y[i];
        }

        fftwf_execute(p);

        for(i=0;i<npts/2+1;i++)
        {
                printf("%f %f\n",i*df,sqrtf(out[i][0]*out[i][0]+out[i][1]*out[i][1]));
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);

        return 0;
}
