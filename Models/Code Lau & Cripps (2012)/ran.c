#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "GARCH.h"

void r_gamma(double *alpha,double *beta,double *r_gamma_s) {
	double *scale;
	scale = (double *) calloc(1,sizeof(double));
	GetRNGstate();
	*scale = 1.0 / (*beta);
	*r_gamma_s = rgamma(*alpha,*scale);
	PutRNGstate();
	free(scale);
	}

void ln_r_unif(double *ln_r_unif_s) {
	GetRNGstate();
	*ln_r_unif_s = log(runif(0.0,1.0));
	PutRNGstate();
	}

void r_unif(double *r_unif_s) {
	GetRNGstate();
	*r_unif_s = runif(0.0,1.0);
	PutRNGstate();
	}

void chol(double *matrixA,int *dimA,double *chol_s) {
	int INFO;
	int i,j;

	for (i=0;i<((*dimA)*(*dimA));i++) { chol_s[i] = matrixA[i]; }
	F77_NAME(dpotrf)("L",dimA,chol_s,dimA,&INFO);
	for (i=0;i<(*dimA);i++) { for (j=0;j<i;j++) { chol_s[i*(*dimA)+j] = 0.0; } }
	if (INFO!=0) { Rprintf("Error in Choleski Decomposition\n"); }
	}

void chol2inv(double *chol_s,int *dimA,double *inv_s) {
	int INFO;
	int i,j;

	for (i=0;i<((*dimA)*(*dimA));i++) { inv_s[i] = chol_s[i]; }
	F77_NAME(dpotri)("L",dimA,inv_s,dimA,&INFO);
	for (i=0;i<(*dimA);i++) { for (j=0;j<i;j++) { inv_s[i*(*dimA)+j] = inv_s[j*(*dimA)+i]; } }
	if (INFO!=0) { Rprintf("Error in Choleski Decomposition to Inverse\n"); }
	}

void inv(double *matrixA,int *dimA,double *inv_s) {
	int INFO;
	int i;

	double *matrixB;
	int *IPIV;
	matrixB = (double *) calloc((*dimA)*(*dimA),sizeof(double));
	IPIV = (int *) calloc((*dimA),sizeof(int));

	for (i=0;i<((*dimA)*(*dimA));i++) { inv_s[i] = matrixA[i]; }
	// Memcpy(invA,matrixA,((*dimA)*(*dimA)));
	for (i=0;i<((*dimA)*(*dimA));i++) { matrixB[i] = 0.0; }
	for (i=0;i<(*dimA);i++) { matrixB[i*(*dimA)+i] = 1.0; }
	F77_NAME(dgesv)(dimA,dimA,inv_s,dimA,IPIV,matrixB,dimA,&INFO);
	for (i=0;i<((*dimA)*(*dimA));i++) { inv_s[i] = matrixB[i]; }
	// Memcpy(matrixB,invA,((*dimA)*(*dimA)));
	if (INFO!=0) { Rprintf("Error in Computing the Inverse\n"); }
	free(matrixB);
	free(IPIV);
	}

void eigen(double *matrixA,int *dimA,double *eigenValuesA,double *eigenVectorA) {
	int INFO;
    int lwork,liwork,itmp,m;
    double *work,tmp;
	int *iwork;
	double vl=0.0,vu=0.0,abstol=0.0;
	int il,iu,i;

	int *isuppz; isuppz = (int *) calloc(2*(*dimA),sizeof(int));
	double *matrixA_s; matrixA_s = (double *) calloc((*dimA)*(*dimA),sizeof(double));
	for (i=0;i<(*dimA)*(*dimA);i++) { matrixA_s[i] = matrixA[i]; }

    lwork = -1;
	liwork = -1;
	F77_NAME(dsyevr)("V","A","L",dimA,matrixA_s,dimA,
		 &vl,&vu,&il,&iu,&abstol,&m,eigenValuesA,
		 eigenVectorA,dimA,isuppz,&tmp,&lwork,&itmp,&liwork,&INFO);

    lwork = (int) tmp;
    liwork = itmp;
    work = (double *) calloc(lwork,sizeof(double));
    iwork = (int *) calloc(liwork,sizeof(int));

	F77_NAME(dsyevr)("V","A","L",dimA,matrixA_s,dimA,
		 &vl,&vu,&il,&iu,&abstol,&m,eigenValuesA,
		 eigenVectorA,dimA,isuppz,work,&lwork,iwork,&liwork,&INFO);

	if (INFO!=0) { Rprintf("Error in Computing the Inverse\n"); }
	free(work);
	free(iwork);
	free(isuppz);
	free(matrixA_s);
	}

void ln_d_t(double *XX,double *meanA,double *sigmaA,double *dfA,double *lnretval_s) {
	double distval = ((*XX)-(*meanA)) * ((*XX)-(*meanA)) / (*sigmaA);
	double logdet = log(*sigmaA);
	*lnretval_s = lgammafn(0.5*((*dfA)+1.0)) - lgammafn(0.5*(*dfA)) 
		- 0.5 * (log((*dfA)*M_PI)+logdet+((*dfA)+1.0)*log(1.0+distval/(*dfA)));
	}

void ln_d_norm(double *XX,double *meanA,double *sigmaA,double *lnretval_s) {
	double distval = ((*XX)-(*meanA)) * ((*XX)-(*meanA)) / (*sigmaA);
	double logdet = log(*sigmaA);
	*lnretval_s = -0.5 * (log(2*PI)+logdet+distval);
	}

void ln_d_uv_norm(double *XX,double *meanA,double *sigmaA,double *lnretval_s) {
	double distval = ((*XX)-(*meanA)) * ((*XX)-(*meanA)) / (*sigmaA);
	double logdet = log(*sigmaA);
	*lnretval_s = -0.5 * (log(2*PI)+logdet+distval);
	}
	
void ln_d_mv_norm(double *XX,double *meanA,double *sigmaA,int *dimA,double *lnretval_s) {
	double distval;
	double logdet;
	int i,j;
	
	double *invA;
	double *tmpvector;
	double *evalues;
	double *evectors;
	invA = (double *) calloc((*dimA)*(*dimA),sizeof(double));
	tmpvector = (double *) calloc((*dimA),sizeof(double));
	evalues = (double *) calloc((*dimA)*(*dimA),sizeof(double));
	evectors = (double *) calloc((*dimA)*(*dimA),sizeof(double));
	
	if (*dimA>1) {
		inv(sigmaA,dimA,invA);
		distval = 0.0;
		for (i=0;i<(*dimA);i++) {
			tmpvector[i] = 0.0;
			for (j=0;j<(*dimA);j++) {
				tmpvector[i] = tmpvector[i] + invA[j*(*dimA)+i] * (XX[j]-meanA[j]);
				}
			distval = distval + tmpvector[i] * (XX[i]-meanA[i]);
			}
		eigen(sigmaA,dimA,evalues,evectors);
		logdet = 0.0;
		for (i=0;i<(*dimA);i++) {
			logdet = logdet + log(evalues[i]);
			}
		} else {
		*tmpvector = (*XX) - (*meanA);
		distval = (*tmpvector) * (*tmpvector) / (*sigmaA);
		logdet = log(*sigmaA);
		}
	*lnretval_s = -0.5 * ((*dimA)*log(2.0*PI)+logdet+distval);

	free(invA);
	free(tmpvector);
	free(evalues);
	free(evectors);
	}

void r_mv_norm_chol(double *meanA,double *cholA,int *dimA,double *r_mv_norm_s) {
	int i,j;
	double *std_norm_s;
	std_norm_s = (double *) calloc((*dimA),sizeof(double));

	GetRNGstate();
	for (i=0;i<(*dimA);i++) { std_norm_s[i] = rnorm(0.0,1.0); }
	for (i=0;i<(*dimA);i++) {
		r_mv_norm_s[i] = meanA[i];
		for (j=0;j<(*dimA);j++) {
			r_mv_norm_s[i] = r_mv_norm_s[i] + cholA[j*(*dimA)+i] * std_norm_s[j];
			}
		}
	PutRNGstate();
	free(std_norm_s);
	}

void r_mv_norm(double *meanA,double *sigmaA,int *dimA,double *r_mv_norm_s) {
	int i,j;
	double *cholA;
	double *std_norm_s;
	cholA = (double *) calloc((*dimA)*(*dimA),sizeof(double));
	std_norm_s = (double *) calloc((*dimA),sizeof(double));

	GetRNGstate();
	chol(sigmaA,dimA,cholA);
	for (i=0;i<(*dimA);i++) { std_norm_s[i] = rnorm(0.0,1.0); }
	// for (i=0;i<(*dimA);i++) { Rprintf("%lf\t",std_norm_s[i]); } Rprintf("\n");
	for (i=0;i<(*dimA);i++) {
		r_mv_norm_s[i] = meanA[i];
		for (j=0;j<(*dimA);j++) {
			r_mv_norm_s[i] = r_mv_norm_s[i] + cholA[j*(*dimA)+i] * std_norm_s[j];
			}
		}
	// for (i=0;i<(*dimA);i++) { for (j=0;j<(*dimA);j++) { Rprintf("%lf\t",cholA[j*(*dimA)+i]); } Rprintf("\n"); }
	PutRNGstate();
	free(cholA);
	free(std_norm_s);
	}

void r_uv_norm(double *meanA,double *sigmaA,int *dimA,double *r_uv_norm_s) {
	GetRNGstate();
	*r_uv_norm_s = *meanA + sqrt(*sigmaA) * rnorm(0.0,1.0);
	PutRNGstate();
	}

void r_norm(double *meanA,double *sigmaA,double *r_uv_norm_s) {
	GetRNGstate();
	*r_uv_norm_s = *meanA + sqrt(*sigmaA) * rnorm(0.0,1.0);
	PutRNGstate();
	}

void r_Ltnorm(double *mean,double *var,double *rtnormA) {
	double sd = sqrt(*var);
	double alpha = -(*mean) / sd;
    if (alpha < 0.45) {
		*rtnormA = (*mean) + sd * nrs_a_inf(&alpha);
		} else {
		*rtnormA = (*mean) + sd * ers_a_inf(&alpha);
		}
	}

double nrs_a_inf(double *a) {
	double x = (*a) - 1.0;
	GetRNGstate();
	while(x<(*a)) { x = rnorm(0.0,1.0); }
	PutRNGstate();
	return x;
	}

double ers_a_inf(double *a) {
	double ainv = 1.0 / (*a);
	double x,rho;
	GetRNGstate();
	do{
		x = rexp(ainv) + (*a); /* rexp works with 1/lambda */
		rho = exp(-0.5*pow((x-(*a)),2));
		} while (runif(0.0,1.0)>rho);
	PutRNGstate();
	return x;
	}

void ln_d_Ltnorm(double *XX,double *mean,double *var,double *lnretval_s) {
	double sd = sqrt(*var);
	double alpha = 0.0;	
	*lnretval_s = dnorm(*XX,*mean,sd,1) - pnorm(alpha,*mean,sd,0,1);
	}

