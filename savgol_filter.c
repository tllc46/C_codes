/*
gcc savgol_filter.c -lm -lgsl -llapack -lrefblas -lgfortran
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct savitzky_golay
{
	int window_length;
	int polyorder;
	int halflen;
	int rem;
	double *coeffs;
	double *left_fit_mat;
	double *right_fit_mat;
} savgol;

void dgelsd_(int *m,int *n,int *nrhs,double *a,int *lda,double *b,int *ldb,double *s,double *rcond,int *rank,double *work,int *lwork,int *iwork,int *info);
void dpotrf_(char *uplo,int *n,double *a,int *lda,int *info);
void dpotri_(char *uplo,int *n,double *a,int *lda,int *info);
void dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
void dgemm_(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
void dgemv_(char *trans,int *m,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);

int Faulhaber(int n,int power)
{
	//important: sum from 0 to n-1
	int sum;

	switch(power)
	{
		case 0:
			sum=n;
			return sum;
		case 1:
			sum=(n-1)*n/2;
			return sum;
		case 2:
			sum=(n-1)*n*(2*n-1)/6;
			return sum;
		case 3:
			sum=(n-1)*(n-1)*n*n/4;
			return sum;
		case 4:
			sum=(n-1)*n*(2*n-1)*(3*n*n-3*n-1)/30;
			return sum;
		case 5:
			sum=(n-1)*(n-1)*n*n*(2*n*n-2*n-1)/12;
			return sum;
		default:
			int i;
			sum=0;
			for(i=1;i<n;i++)
			{
				sum+=pow(i,power);
			}
			return sum;
	}
}

void term_savgol(savgol *x_savgol)
{
	free(x_savgol->coeffs);
	free(x_savgol->left_fit_mat);
	free(x_savgol->right_fit_mat);
}

void init_polyfit(savgol *x_savgol)
{
	//solve A*x=b -> (A^t*A)*x=A^t*b -> x=(A^t*A)^(-1)*A^t*b
	char uplo='U';
	int n,lda;
	double *a;
	int info;
	char side='L';
	int m,ldb,ldc;
	double alpha=1,beta=0;
	double *b,*c;
	char transa='N',transb='N';
	int k;
	double *a_edge;

	int i,j;

	//1. calulate (A^t*A)^(-1)
	n=x_savgol->polyorder+1;
	lda=n;
	a=(double *)malloc(lda*n*sizeof(double)); //A^t*A, (polyorder+1) x (polyorder+1)
	for(j=0;j<n;j++)
	{
		for(i=0;i<=j;i++)
		{
			a[j*lda+i]=Faulhaber(x_savgol->window_length,i+j); //column-wise order
		}
	}

	dpotrf_(&uplo,&n,a,&lda,&info); //Cholesky factorization
	dpotri_(&uplo,&n,a,&lda,&info); //inversion using triangular matrix

	//2. calculate (A^t*A)^(-1)*A^t
	m=x_savgol->polyorder+1;
	n=x_savgol->window_length;
	lda=m;
	ldb=m;
	ldc=m;
	b=(double *)malloc(ldb*n*sizeof(double)); //A^t, (polyorder+1) x window_length
	for(j=0;j<n;j++)
	{
		for(i=0;i<ldb;i++)
		{
			b[j*ldb+i]=pow(j,i); //column-wise order
		}
	}
	c=(double *)malloc(ldc*n*sizeof(double)); //(A^t*A)^(-1)*A^t, (polyorder+1) x window_length

	dsymm_(&side,&uplo,&m,&n,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);

	//3. calculate edge*(A^t*A)^(-1)*A^t
	m=x_savgol->halflen;
	k=x_savgol->polyorder+1;
	n=x_savgol->window_length;
	lda=m;
	a_edge=(double *)malloc(lda*k*sizeof(double)); //halflen x (polyorder+1)
	ldb=k;
	ldc=m;
	x_savgol->left_fit_mat=(double *)malloc(ldc*n*sizeof(double)); //halflen x window_length
	x_savgol->right_fit_mat=(double *)malloc(ldc*n*sizeof(double));

	//left edge
	for(j=0;j<k;j++)
	{
		for(i=0;i<lda;i++)
		{
			a_edge[j*lda+i]=pow(i,j); //column-wise order
		}
	}
	dgemm_(&transa,&transb,&m,&n,&k,&alpha,a_edge,&lda,c,&ldb,&beta,x_savgol->left_fit_mat,&ldc);

	//right edge
	for(j=0;j<k;j++)
	{
		for(i=0;i<lda;i++)
		{
			a_edge[j*lda+i]=pow((n-x_savgol->halflen)+i,j); //column-wise order
		}
	}
	dgemm_(&transa,&transb,&m,&n,&k,&alpha,a_edge,&lda,c,&ldb,&beta,x_savgol->right_fit_mat,&ldc);

	free(a);
	free(b);
	free(c);
	free(a_edge);
}

void init_savgol(savgol *x_savgol,int window_length,int polyorder)
{
	double pos,x;
	int m,n,nrhs=1,lda,ldb;
	double *A,*s;
	double rcond=-1;
	double *work=(double *)malloc(1*sizeof(double));
	int rank,lwork=-1,info,*iwork=(int *)malloc(1*sizeof(int));

	int i,j;

	//initialize savgol structure
	x_savgol->window_length=window_length;
	x_savgol->polyorder=polyorder;
	x_savgol->halflen=x_savgol->window_length/2;
	x_savgol->rem=x_savgol->window_length%2;

	if(x_savgol->rem==0)
	{
		pos=x_savgol->halflen-0.5;
	}
	else
	{
		pos=x_savgol->halflen;
	}

	m=x_savgol->polyorder+1;
	n=x_savgol->window_length;
	lda=m;
	ldb=n;
	A=(double *)malloc(lda*n*sizeof(double)); //(polyorder+1) x window_length
	for(j=0;j<n;j++)
	{
		x=pos-j;
		for(i=0;i<lda;i++)
		{
			A[j*lda+i]=pow(x,i); //columnwise order
		}
	}
	x_savgol->coeffs=(double *)malloc(ldb*sizeof(double)); //window_length x 1
	memset(x_savgol->coeffs,0,m*sizeof(int));
	x_savgol->coeffs[0]=1;
	s=(double *)malloc(m*sizeof(double));

	dgelsd_(&m,&n,&nrhs,A,&lda,x_savgol->coeffs,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info); //workspace query

	lwork=work[0];
	work=(double *)malloc((int)lwork*sizeof(double));
	iwork=(int *)malloc(iwork[0]*sizeof(int));

	dgelsd_(&m,&n,&nrhs,A,&lda,x_savgol->coeffs,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info); //main calculation

	init_polyfit(x_savgol);

	free(A);
	free(s);
	free(work);
	free(iwork);
}

void convolve1d(savgol *x_savgol,double *data,int npts_data,double *data_cnvlv)
{
	int start_idx;

	int i,j;

	memset(data_cnvlv+x_savgol->halflen,0,(npts_data-2*x_savgol->halflen)*sizeof(double));

	//convolve data and coeffs
	for(i=x_savgol->halflen;i<npts_data-x_savgol->halflen;i++)
	{
		start_idx=i-(x_savgol->window_length-1)/2;
		for(j=0;j<x_savgol->window_length;j++)
		{
			data_cnvlv[i]+=data[start_idx+j]*x_savgol->coeffs[j];
		}
	}
}

void savgol_filter(savgol *x_savgol,double *data,int npts_data,double *data_fltr)
{
	char trans='N';
	int m=x_savgol->halflen,n=x_savgol->window_length,lda=m,incx=1,incy=1;
	double alpha=1,beta=0;
	double *data_fit=(double *)malloc(m*sizeof(double));

	convolve1d(x_savgol,data,npts_data,data_fltr);

	//left edge fit
	dgemv_(&trans,&m,&n,&alpha,x_savgol->left_fit_mat,&lda,data,&incx,&beta,data_fit,&incy);
	memmove(data_fltr,data_fit,m*sizeof(double));

	//right edge fit
	dgemv_(&trans,&m,&n,&alpha,x_savgol->right_fit_mat,&lda,data+(npts_data-n),&incx,&beta,data_fit,&incy);
	memmove(data_fltr+(npts_data-m),data_fit,m*sizeof(double));
}

int main(int argc,char **argv)
{
	savgol x_savgol;
	double data[21]={1,4,3,5,8,3,2,4,5,8,9,3,5,3,2,1,4,2,6,8,3};
	double *data_fltr=(double *)malloc(21*sizeof(double));

	int i;

	init_savgol(&x_savgol,10,2);
	savgol_filter(&x_savgol,data,21,data_fltr);
	for(i=0;i<21;i++)
	{
		printf("data_fltr[%d]=%f\n",i,data_fltr[i]);
	}

	term_savgol(&x_savgol);

	return 0;
}
