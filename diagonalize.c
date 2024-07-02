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
        int n=3,lda=3,lwork,info,i,j;
        double *w,*rwork;
        f_complex *work;
        f_complex a[3][3]={ //ForTran에 전달해야 하기 때문에 열 방향 정렬
                {{6,0},{2,-8},{9,5}},
                {{2,8},{4,0},{7,3}},
                {{9,-5},{7,-3},{1,0}}
        };
        /* 실제 행렬
        6 2+8i 9-5i
        2-8i 4 7-3i
        9+5i 7+3i 1
        */
        char jobz='V',uplo='U';

        w=(double *)malloc(sizeof(double)*n);
        rwork=(double *)malloc(sizeof(double)*(3*n-2));
        lwork=-1;
        work=(f_complex *)malloc(sizeof(f_complex)); //최적 lwork만 담기 위해 크기 1

        zheev_(&jobz,&uplo,&n,a,&lda,w,work,&lwork,rwork,&info); //최적 lwork 추론

        lwork=work[0].real; //최적 lwork
        work=(f_complex *)malloc(sizeof(f_complex)*lwork); //다시 크기 lwork만큼 할당

        zheev_(&jobz,&uplo,&n,a,&lda,w,work,&lwork,rwork,&info); //실제 계산

        printf("info=%d\n\n",info);
        printf("eigenvalues");
        for(i=0;i<n;i++)
        {
                printf("w[%d]=%f\n",i,w[i]); //오름차순 고윳값
        }
        printf("\n");
        for(i=0;i<n;i++)
        {
                printf("eigenvector %d\n",i);
                for(j=0;j<n;j++)
                {
                        printf("a[%d][%d]=%.9f + %.9fi\n",j,i,a[i][j].real,a[i][j].imag);
                }
        }

        free(w);
        free(rwork);
        free(work);

        return 0;
}
