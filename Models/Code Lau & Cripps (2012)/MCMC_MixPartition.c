#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "GARCH.h"

void sample_PD_GIB_1STEP_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *X,			// [7]
	double *SUM_X,		// [8]
	double *SOS_X		// [9]
	) { 

	int i,j,index_c;
	// int nn,kk;
	// double lnV1,lnV2;
	double ln_r_unif_s;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *PD_alpha; PD_alpha = (double *) calloc(1,sizeof(double));
	double *PD_q; PD_q = (double *) calloc(1,sizeof(double));
	double p_a,p_b,p_m,p_t,ln_d_t_s,mean_t,vari_t,df_t;

	*PD_alpha = *PARAs;
	*PD_q = *(PARAs+3);

	for (i=0;i<*Nsize;i++) {

		index_c = C[i] - 1;
		e[index_c] = e[index_c] - 1;
		SUM_X[index_c] = SUM_X[index_c] - X[i];
		SOS_X[index_c] = SOS_X[index_c] - (X[i]*X[i]);
		if (e[index_c]<1) {
			for (j=0;j<(*Nsize);j++) { if (C[j]==(*np)) { C[j] = C[i]; } }
			*np = *np - 1;
			e[index_c] = e[*np];
			SUM_X[index_c] = SUM_X[*np];
			SOS_X[index_c] = SOS_X[*np];
			}
		e[*np] = 0;
		SUM_X[*np] = 0.0;
		SOS_X[*np] = 0.0;

		for (j=0;j<(*np);j++) {
			p_t = PRIORs[3] + e[j];
			p_m = (PRIORs[3]*PRIORs[2]+SUM_X[j]) / p_t;
			p_a = PRIORs[0] + 0.5 * e[j];
			p_b = PRIORs[1] + 0.5 * (SOS_X[j]-SUM_X[j]*SUM_X[j]/e[j])
				+ 0.5 * (PRIORs[2]-SUM_X[j]/e[j]) * (PRIORs[2]-SUM_X[j]/e[j])
				/ (1.0/PRIORs[3]+ 1.0/e[j]);
			mean_t = p_m;
			vari_t = (p_t+1.0) / p_t * p_b / p_a;
			df_t = 2.0 * p_a;
			ln_d_t(X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
			lnW[j] = log(e[j]-(*PD_alpha)) + ln_d_t_s;
			}
		p_a = PRIORs[0];
		p_b = PRIORs[1];
		p_m = PRIORs[2];
		p_t = PRIORs[3];
		mean_t = p_m;
		vari_t = (p_t+1.0) / p_t * p_b / p_a;
		df_t = 2.0 * p_a;
		ln_d_t(X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
		lnW[*np] = log((*PD_q)+(*PD_alpha)*(*np)) + ln_d_t_s;

		for (j=1;j<(*np+1);j++) {
			if (lnW[j-1]>lnW[j]) {
				lnW[j] = lnW[j-1] + log(1.0+exp(lnW[j]-lnW[j-1]));
				} else {
				lnW[j] = lnW[j] + log(1.0+exp(lnW[j-1]-lnW[j]));
				}
			}
		for (j=0;j<(*np+1);j++) { lnW[j] = lnW[j] - lnW[*np]; }

		// for (j=0;j<(*np+1);j++) { Rprintf("%lf\t",lnW[j]); } Rprintf("\n");
		
		ln_r_unif(&ln_r_unif_s);
		index_c = *np;
		for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { index_c = j; j = *np; } }
		C[i] = index_c + 1;
		e[index_c] = e[index_c] + 1;
		SUM_X[index_c] = SUM_X[index_c] + X[i];
		SOS_X[index_c] = SOS_X[index_c] + (X[i]*X[i]);
		if (index_c==(*np)) { *np = *np + 1; }
		// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",C[j]); } Rprintf("\n");
		// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",e[j]); } Rprintf("\n");
		// for (j=0;j<*Nsize;j++) { Rprintf("%lf\t",lnW[j]); } Rprintf("\n");
		}
	free(PD_alpha);
	free(PD_q);
	free(lnW);

	}

void sample_PD_GIB_2STEP_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *X,			// [7]
	double *SUM_X,		// [8]
	double *SOS_X		// [9]
	) { 

	int i,j,index_c;
	double ln_r_unif_s;
	double gamma_a,gamma_b,U_alpha;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *PD_alpha; PD_alpha = (double *) calloc(1,sizeof(double)); *PD_alpha = *PARAs;
	double *PD_q; PD_q = (double *) calloc(1,sizeof(double)); *PD_q = *(PARAs+3);
	double p_a,p_b,p_m,p_t,ln_d_t_s,mean_t,vari_t,df_t;

	for (i=0;i<*Nsize;i++) {

		gamma_a = *PD_q / (*PD_alpha) + (*np);
		gamma_b = 1 / (*PD_alpha);
		r_gamma(&gamma_a,&gamma_b,&U_alpha); // U = exp(log(U_alpha)/(*PD_alpha));
		
		index_c = C[i] - 1;
		e[index_c] = e[index_c] - 1;
		SUM_X[index_c] = SUM_X[index_c] - X[i];
		SOS_X[index_c] = SOS_X[index_c] - (X[i]*X[i]);
		if (e[index_c]<1) {
			for (j=0;j<(*Nsize);j++) { if (C[j]==(*np)) { C[j] = C[i]; } }
			*np = *np - 1;
			e[index_c] = e[*np];
			SUM_X[index_c] = SUM_X[*np];
			SOS_X[index_c] = SOS_X[*np];
			}
		e[*np] = 0;
		SUM_X[*np] = 0.0;
		SOS_X[*np] = 0.0;

		for (j=0;j<(*np);j++) {
			p_t = PRIORs[3] + e[j];
			p_m = (PRIORs[3]*PRIORs[2]+SUM_X[j]) / p_t;
			p_a = PRIORs[0] + 0.5 * e[j];
			p_b = PRIORs[1] + 0.5 * (SOS_X[j]-SUM_X[j]*SUM_X[j]/e[j])
				+ 0.5 * (PRIORs[2]-SUM_X[j]/e[j]) * (PRIORs[2]-SUM_X[j]/e[j])
				/ (1.0/PRIORs[3]+ 1.0/e[j]);
			mean_t = p_m;
			vari_t = (p_t+1.0) / p_t * p_b / p_a;
			df_t = 2.0 * p_a;
			ln_d_t(X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
			lnW[j] = log(e[j]-(*PD_alpha)) + ln_d_t_s;
			}
		p_a = PRIORs[0];
		p_b = PRIORs[1];
		p_m = PRIORs[2];
		p_t = PRIORs[3];
		mean_t = p_m;
		vari_t = (p_t+1.0) / p_t * p_b / p_a;
		df_t = 2.0 * p_a;
		ln_d_t(X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
		lnW[*np] = log(U_alpha) + ln_d_t_s;

		for (j=1;j<(*np+1);j++) {
			if (lnW[j-1]>lnW[j]) {
				lnW[j] = lnW[j-1] + log(1.0+exp(lnW[j]-lnW[j-1]));
				} else {
				lnW[j] = lnW[j] + log(1.0+exp(lnW[j-1]-lnW[j]));
				}
			}
		for (j=0;j<(*np+1);j++) { lnW[j] = lnW[j] - lnW[*np]; }

		ln_r_unif(&ln_r_unif_s);
		index_c = *np;
		for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { index_c = j; j = *np; } }
		C[i] = index_c + 1;
		e[index_c] = e[index_c] + 1;
		SUM_X[index_c] = SUM_X[index_c] + X[i];
		SOS_X[index_c] = SOS_X[index_c] + (X[i]*X[i]);
		if (index_c==(*np)) { *np = *np + 1; }
		// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",C[j]); } Rprintf("\n");
		// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",e[j]); } Rprintf("\n");
		// for (j=0;j<*Nsize;j++) { Rprintf("%lf\t",lnW[j]); } Rprintf("\n");
		}
	free(PD_alpha);
	free(PD_q);
	free(lnW);
	}

void sample_NGG_GIB_2STEP_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *X,			// [7]
	double *SUM_X,		// [8]
	double *SOS_X		// [9]
	) {

	int i,j,index_c;
	double ln_r_unif_s;
	double gamma_a,gamma_b,r_gamma_s,logpsi,V;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *NGG_alpha; NGG_alpha = (double *) calloc(1,sizeof(double));
	double *NGG_theta; NGG_theta = (double *) calloc(1,sizeof(double));
	double *NGG_b; NGG_b = (double *) calloc(1,sizeof(double));
	double p_a,p_b,p_m,p_t,ln_d_t_s,mean_t,vari_t,df_t;
	
	*NGG_alpha = *PARAs;
	*NGG_theta = *(PARAs+1);
	*NGG_b = *(PARAs+2);

	for (i=0;i<*Nsize;i++) {

		// std::cout << " in " << std::endl;
		do {
			gamma_a = *np;
			gamma_b = *NGG_theta / *NGG_alpha;
			r_gamma(&gamma_a,&gamma_b,&r_gamma_s);
			V = exp(log(r_gamma_s)/(*NGG_alpha));
			// V = pow(r_gamma_s,1.0/(*NGG_alpha));
			logpsi = -((*NGG_theta)/(*NGG_alpha)) * (pow(V+*NGG_b,*NGG_alpha)-pow(V,*NGG_alpha))
				+ (*Nsize-(*NGG_alpha)*(*np))*(log(V)-log(V+*NGG_b));
			// logpsi = -((*NGG_theta)/(*NGG_alpha)*(exp((*NGG_alpha)*log(V+*NGG_b))-exp((*NGG_alpha)*log(V))))
			// 	+ (*Nsize-(*NGG_alpha)*(*np)+1.0)*(log(V)-log(V+*NGG_b));
			ln_r_unif(&ln_r_unif_s);
			// std::cout << log(runif) << ' ' << logpsi << std::endl;
			} while (ln_r_unif_s>logpsi); // U = V;
		// std::cout << " out " << std::endl;

		index_c = C[i] - 1;
		e[index_c] = e[index_c] - 1;
		SUM_X[index_c] = SUM_X[index_c] - X[i];
		SOS_X[index_c] = SOS_X[index_c] - (X[i]*X[i]);
		if (e[index_c]<1) {
			for (j=0;j<(*Nsize);j++) { if (C[j]==(*np)) { C[j] = C[i]; } }
			*np = *np - 1;
			e[index_c] = e[*np];
			SUM_X[index_c] = SUM_X[*np];
			SOS_X[index_c] = SOS_X[*np];
			}
		e[*np] = 0;
		SUM_X[*np] = 0.0;
		SOS_X[*np] = 0.0;

		for (j=0;j<(*np);j++) {
			p_t = PRIORs[3] + e[j];
			p_m = (PRIORs[3]*PRIORs[2]+SUM_X[j]) / p_t;
			p_a = PRIORs[0] + 0.5 * e[j];
			p_b = PRIORs[1] + 0.5 * (SOS_X[j]-SUM_X[j]*SUM_X[j]/e[j])
				+ 0.5 * (PRIORs[2]-SUM_X[j]/e[j]) * (PRIORs[2]-SUM_X[j]/e[j])
				/ (1.0/PRIORs[3]+ 1.0/e[j]);
			mean_t = p_m;
			vari_t = (p_t+1.0) / p_t * p_b / p_a;
			df_t = 2.0 * p_a;
			ln_d_t(X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
			lnW[j] = log(e[j]-(*NGG_alpha)) + ln_d_t_s;
			}
		p_a = PRIORs[0];
		p_b = PRIORs[1];
		p_m = PRIORs[2];
		p_t = PRIORs[3];
		mean_t = p_m;
		vari_t = (p_t+1.0) / p_t * p_b / p_a;
		df_t = 2.0 * p_a;
		ln_d_t(X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
		lnW[*np] = log(*NGG_theta) + (*NGG_alpha) * log(V+(*NGG_b)) + ln_d_t_s;

		for (j=1;j<(*np+1);j++) {
			if (lnW[j-1]>lnW[j]) {
				lnW[j] = lnW[j-1] + log(1.0+exp(lnW[j]-lnW[j-1]));
				} else {
				lnW[j] = lnW[j] + log(1.0+exp(lnW[j-1]-lnW[j]));
				}
			}
		for (j=0;j<(*np+1);j++) { lnW[j] = lnW[j] - lnW[*np]; }

		ln_r_unif(&ln_r_unif_s);
		index_c = *np;
		for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { index_c = j; j = *np; } }
		C[i] = index_c + 1;
		e[index_c] = e[index_c] + 1;
		SUM_X[index_c] = SUM_X[index_c] + X[i];
		SOS_X[index_c] = SOS_X[index_c] + (X[i]*X[i]);
		if (index_c==(*np)) { *np = *np + 1; }
		// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",C[j]); } Rprintf("\n");
		// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",e[j]); } Rprintf("\n");
		// for (j=0;j<*Nsize;j++) { Rprintf("%lf\t",lnW[j]); } Rprintf("\n");
		}
	free(NGG_alpha);
	free(NGG_theta);
	free(NGG_b);
	free(lnW);
	}

void initial_GIB_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *X,			// [7]
	double *SOS_X,		// [8]
	double *SUM_X		// [9]
	) {
	int i;
	for (i=0;i<(*Nsize);i++) {
		C[i] = i + 1;
		e[i] = 1;
		SUM_X[i] = X[i];
		SOS_X[i] = X[i] * X[i];
		}
	*np = *Nsize;
	// for (i=0;i<*Nsize;i++) { Rprintf("%d\t",C[i]); } Rprintf("\n");
	// for (i=0;i<*Nsize;i++) { Rprintf("%d\t",e[i]); } Rprintf("\n");
	}

void Pred_PD_GIB_1STEP_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *SUM_X,		// [7]
	double *SOS_X,		// [8]
	double *pred_X,		// [9]
	double *pred_lnY,	// [10]
	int *pred_n			// [11]
	) { 

	int i,j;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *PD_alpha; PD_alpha = (double *) calloc(1,sizeof(double));
	double *PD_q; PD_q = (double *) calloc(1,sizeof(double));
	double p_a,p_b,p_m,p_t,ln_d_t_s,mean_t,vari_t,df_t;

	*PD_alpha = *PARAs;
	*PD_q = *(PARAs+3);

	for (i=0;i<*pred_n;i++) {

		for (j=0;j<(*np);j++) {
			p_t = PRIORs[3] + e[j];
			p_m = (PRIORs[3]*PRIORs[2]+SUM_X[j]) / p_t;
			p_a = PRIORs[0] + 0.5 * e[j];
			p_b = PRIORs[1] + 0.5 * (SOS_X[j]-SUM_X[j]*SUM_X[j]/e[j])
				+ 0.5 * (PRIORs[2]-SUM_X[j]/e[j]) * (PRIORs[2]-SUM_X[j]/e[j])
				/ (1.0/PRIORs[3]+ 1.0/e[j]);
			mean_t = p_m;
			vari_t = (p_t+1.0) / p_t * p_b / p_a;
			df_t = 2.0 * p_a;
			ln_d_t(pred_X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
			lnW[j] = log(e[j]-(*PD_alpha)) - log(*Nsize+(*PD_q)) + ln_d_t_s;
			}
		p_a = PRIORs[0];
		p_b = PRIORs[1];
		p_m = PRIORs[2];
		p_t = PRIORs[3];
		mean_t = p_m;
		vari_t = (p_t+1.0) / p_t * p_b / p_a;
		df_t = 2.0 * p_a;
		ln_d_t(pred_X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
		lnW[*np] = log((*PD_q)+(*PD_alpha)*(*np)) - log(*Nsize+(*PD_q)) + ln_d_t_s;

		for (j=1;j<(*np+1);j++) {
			if (lnW[j-1]>lnW[j]) {
				lnW[j] = lnW[j-1] + log(1.0+exp(lnW[j]-lnW[j-1]));
				} else {
				lnW[j] = lnW[j] + log(1.0+exp(lnW[j-1]-lnW[j]));
				}
			}

		pred_lnY[i] = lnW[*np];
		}
	}

void Pred_PD_GIB_2STEP_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *SUM_X,		// [7]
	double *SOS_X,		// [8]
	double *pred_X,		// [9]
	double *pred_lnY,	// [10]
	int *pred_n			// [11]
	) { 

	int i,j;
	double gamma_a,gamma_b,U_alpha;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *PD_alpha; PD_alpha = (double *) calloc(1,sizeof(double));
	double *PD_q; PD_q = (double *) calloc(1,sizeof(double));
	double p_a,p_b,p_m,p_t,ln_d_t_s,mean_t,vari_t,df_t;

	*PD_alpha = *PARAs;
	*PD_q = *(PARAs+3);

	for (i=0;i<*pred_n;i++) {

		gamma_a = *PD_q / (*PD_alpha) + (*np);
		gamma_b = 1 / (*PD_alpha);
		r_gamma(&gamma_a,&gamma_b,&U_alpha);

		for (j=0;j<(*np);j++) {
			p_t = PRIORs[3] + e[j];
			p_m = (PRIORs[3]*PRIORs[2]+SUM_X[j]) / p_t;
			p_a = PRIORs[0] + 0.5 * e[j];
			p_b = PRIORs[1] + 0.5 * (SOS_X[j]-SUM_X[j]*SUM_X[j]/e[j])
				+ 0.5 * (PRIORs[2]-SUM_X[j]/e[j]) * (PRIORs[2]-SUM_X[j]/e[j])
				/ (1.0/PRIORs[3]+ 1.0/e[j]);
			mean_t = p_m;
			vari_t = (p_t+1.0) / p_t * p_b / p_a;
			df_t = 2.0 * p_a;
			ln_d_t(pred_X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
			lnW[j] = log(e[j]-(*PD_alpha)) - log((*Nsize)-(*np)*(*PD_alpha)+U_alpha) + ln_d_t_s;
			}
		p_a = PRIORs[0];
		p_b = PRIORs[1];
		p_m = PRIORs[2];
		p_t = PRIORs[3];
		mean_t = p_m;
		vari_t = (p_t+1.0) / p_t * p_b / p_a;
		df_t = 2.0 * p_a;
		ln_d_t(pred_X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
		lnW[*np] = log(U_alpha) - log((*Nsize)-(*np)*(*PD_alpha)+U_alpha) + ln_d_t_s;

		for (j=1;j<(*np+1);j++) {
			if (lnW[j-1]>lnW[j]) {
				lnW[j] = lnW[j-1] + log(1.0+exp(lnW[j]-lnW[j-1]));
				} else {
				lnW[j] = lnW[j] + log(1.0+exp(lnW[j-1]-lnW[j]));
				}
			}

		pred_lnY[i] = lnW[*np];
		}
	}

void Pred_NGG_GIB_2STEP_MN(
	int *C, 			// [1]
	int *e, 			// [2]
	int *Nsize, 		// [3]
	int *np, 			// [4]
	double *PARAs,		// [5]
	double *PRIORs,		// [6]
	double *SUM_X,		// [7]
	double *SOS_X,		// [8]
	double *pred_X,		// [9]
	double *pred_lnY,	// [10]
	int *pred_n			// [11]
	) { 

	int i,j;
	double gamma_a,gamma_b,r_gamma_s,logpsi,V;
	double new_lnW;
	double ln_r_unif_s;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *NGG_alpha; NGG_alpha = (double *) calloc(1,sizeof(double));
	double *NGG_theta; NGG_theta = (double *) calloc(1,sizeof(double));
	double *NGG_b; NGG_b = (double *) calloc(1,sizeof(double));
	double p_a,p_b,p_m,p_t,ln_d_t_s,mean_t,vari_t,df_t;

	*NGG_alpha = *PARAs;
	*NGG_theta = *(PARAs+1);
	*NGG_b = *(PARAs+2);

	for (i=0;i<*pred_n;i++) {

		do {
			gamma_a = *np;
			gamma_b = *NGG_theta / *NGG_alpha;
			r_gamma(&gamma_a,&gamma_b,&r_gamma_s);
			V = exp(log(r_gamma_s)/(*NGG_alpha));
			logpsi = -((*NGG_theta)/(*NGG_alpha)) * (pow(V+*NGG_b,*NGG_alpha)-pow(V,*NGG_alpha))
				+ (*Nsize-(*NGG_alpha)*(*np))*(log(V)-log(V+*NGG_b));
			ln_r_unif(&ln_r_unif_s);
			} while (ln_r_unif_s>logpsi);

		new_lnW = log(*NGG_theta) + (*NGG_alpha) * log(V+(*NGG_b));

		for (j=0;j<(*np);j++) {
			p_t = PRIORs[3] + e[j];
			p_m = (PRIORs[3]*PRIORs[2]+SUM_X[j]) / p_t;
			p_a = PRIORs[0] + 0.5 * e[j];
			p_b = PRIORs[1] + 0.5 * (SOS_X[j]-SUM_X[j]*SUM_X[j]/e[j])
				+ 0.5 * (PRIORs[2]-SUM_X[j]/e[j]) * (PRIORs[2]-SUM_X[j]/e[j])
				/ (1.0/PRIORs[3]+ 1.0/e[j]);
			mean_t = p_m;
			vari_t = (p_t+1.0) / p_t * p_b / p_a;
			df_t = 2.0 * p_a;
			ln_d_t(pred_X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
			lnW[j] = log(e[j]-(*NGG_alpha)) - log((*Nsize)-(*np)*(*NGG_alpha)+exp(new_lnW)) + ln_d_t_s;
			}
		p_a = PRIORs[0];
		p_b = PRIORs[1];
		p_m = PRIORs[2];
		p_t = PRIORs[3];
		mean_t = p_m;
		vari_t = (p_t+1.0) / p_t * p_b / p_a;
		df_t = 2.0 * p_a;
		ln_d_t(pred_X+i,&mean_t,&vari_t,&df_t,&ln_d_t_s);
		lnW[*np] = new_lnW - log((*Nsize)-(*np)*(*NGG_alpha)+exp(new_lnW)) + ln_d_t_s;

		for (j=1;j<(*np+1);j++) {
			if (lnW[j-1]>lnW[j]) {
				lnW[j] = lnW[j-1] + log(1.0+exp(lnW[j]-lnW[j-1]));
				} else {
				lnW[j] = lnW[j] + log(1.0+exp(lnW[j-1]-lnW[j]));
				}
			}

		pred_lnY[i] = lnW[*np];
		}
	free(lnW);
	free(NGG_alpha);
	free(NGG_theta);
	free(NGG_b);
	}

