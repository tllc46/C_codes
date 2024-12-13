/*
gcc sacpz.c -I/usr/local/src/sac-102.0/inc $(bash sac-config -c -l sac sacio) -lm
*/

#include <stdio.h>
#include <string.h>

#include <sac.h>
#include <sacio.h>
#include <libpz.h> // (source)/inc/libpz.h

pz_t *polezero_try_read(char *pzfile,char *id,char *when); // (source)/icm/libpz.c
Response *response_from_polezero(pz_t *pz,Transfer *t); // (source)/icm/libpz.c

int main(int argc,char **argv)
{
        int max=300000;
        int npts,nerr;
        float y[300000],b,delta;
        char filer[25]="KS.CGYA.HGE.2019.005.SAC";
        int i;

        pz_t *pz; // (binary)/include/sac.h, (source)/inc/icm.h
        Transfer trans; // (source)/inc/libpz.h
        Response *resp; // (source)/inc/libpz.h

        rsac1(filer,y,&npts,&b,&delta,&max,&nerr,strlen(filer));

        pz=polezero_try_read("sacpz_file","KS.CGYA..HGE","2019-01-05"); // (source)/icm/libpz.c

        printf("nzero=%d\n",pz->nzero);
        printf("npole=%d\n",pz->npole);
        for(i=0;i<pz->nzero;i++)
        {
                printf("zero %d=%f %f\n",i+1,pz->zeros[i].re,pz->zeros[i].im);
        }
        for(i=0;i<pz->npole;i++)
        {
                printf("pole %d=%f %f\n",i+1,pz->poles[i].re,pz->poles[i].im);
        }
        printf("constant=%f\n",pz->constant);
        printf("line=%s\n",pz->line);
        printf("\n");

        trans.n=36000;
        trans.dt=delta;
        trans.nfft=trans.n;
        trans.nfreqs=trans.nfft/2+1;
        trans.df=dfreq(trans.nfft,trans.dt);
        trans.f=fftfreq(trans.nfft,0,trans.df);

        printf("n=%d\n",trans.n);
        printf("dt=%f\n",trans.dt);
        printf("nfft=%d\n",trans.nfft);
        printf("nfreqs=%d\n",trans.nfreqs);
        printf("df=%f\n",trans.df);
        printf("\n");

        resp=response_from_polezero(pz,&trans); // (source)/icm/libpz.c

        printf("n=%d\n",resp->n);
        for(i=0;i<5;i++)
        {
                printf("%f + %f\n",resp->re[i],resp->im[i]);
        }

        return 0;
}
