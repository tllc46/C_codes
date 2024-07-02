/* 컴파일 방법
gcc diagonalize.c lapack/liblapack.a lapack/librefblas.a -lgfortran -lm

(1)2차원 배열을 배열 포인터로 전달하기
void zheev_(...,f_complex (*a)[],...);
f_complex a[3][3]=...;
zheev_(...,a,...);

(2)2차원 배열을 포인터로 전달하기
void zheev_(...,f_complex *a,...);
f_complex a[3][3]=...;
zheev_(...,*a,...); 또는 zheev_(...,(f_complex *)a,...);

(3)1차원 배열을 포인터로 전달하기
void zheev_(...,f_complex *a,...);
f_complex a[3*3]=...;
zheev_(...,a,...);

ForTran의 함수, subroutine은 다차원 배열이라도 무조건 첫 번째 원소를 가리키는 포인터로만 받기 때문에
배열 포인터를 전달하더라도 포인터로 인식
*/

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
        double real;
        double imag;
} f_complex;

void zheev_(char *jobz,char *uplo,int *n,f_complex (*a)[],int *lda,double *w,f_complex *work,int *lwork,double *rwork,int *info);

int main(int argc,char **argv)
{
        int n=3,lda=3,lwork=-1,info,i,j;
        double *w,*rwork;
        f_complex *work;
        f_complex a[3][3]={
                {{6,0},{2,-8},{9,5}},
                {{2,8},{4,0},{7,3}},
                {{9,-5},{7,-3},{1,0}}
        };
        char jobz='N',uplo='U';

        w=(double *)malloc(sizeof(double)*n);
        rwork=(double *)malloc(sizeof(double)*(3*n-2));
        work=(f_complex *)malloc(sizeof(f_complex)*5);

        zheev_(&jobz,&uplo,&n,a,&lda,w,work,&lwork,rwork,&info);
        printf("jobz=%c\n",jobz);
        printf("uplo=%c\n",uplo);
        printf("n=%d\n",n);
        printf("lda=%d\n",lda);
        printf("lwork=%d\n",lwork);
        printf("info=%d\n",info);
        for(i=0;i<n;i++)
        {
                printf("w[%d]=%f\n",i,w[i]);
        }
        for(i=0;i<3*n-2;i++)
        {
                printf("rwork[%d]=%f\n",i,rwork[i]);
        }
        for(i=0;i<5;i++)
        {
                printf("work[%d]=%f + %fi\n",i,work[i].real,work[i].imag);
        }
        for(i=0;i<n;i++)
        {
                for(j=0;j<n;j++)
                {
                        printf("a[%d][%d]=%f + %fi\n",i,j,a[j][i].real,a[j][i].imag);
                }
        }

        free(w);
        free(rwork);
        free(work);

        return 0;
}