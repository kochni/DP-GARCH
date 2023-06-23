#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "GARCH.h"

void logLikelihood_i_to_n(
	int *C, 				// [1]
	int *e, 				// [2]
	int *Nsize, 			// [3]
	int *np, 				// [4]
	int *orders,			// [5]  // (max.order,ar.order,garch.a.order,garch.b.order)
	int *dims,				// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARA_ARs,		// [7]  // rep(c(0,rep(0,ar.order)),n.size)
	double *PARA_GARCHs,	// [8]  // rep(c(1,rep(0,garch.a.order),rep(0,garch.b.order)),n.size)
	double *Y_series,		// [9]  // n.size + max.order
	double *m_series,		// [10] // n.size + max.order
	double *e_series,		// [11] // n.size + max.order
	double *h2_series,		// [12] // n.size + max.order
	double *logLikelihood_s,// [13]
	int *ii					// [14]
	) {
	int i,j,index_c;
	int order_max = orders[0];
	int dim_AR1 = dims[0];
	int dim_GARCH1 = dims[2];
	int dim_GARCH_a1 = dims[4];
	int dim_GARCH_b1 = dims[6];
	double ln_d_norm_s;

	*logLikelihood_s = 0.0;
	for (i=*ii;i<*Nsize;i++) {
		index_c = C[i] - 1;
		m_series[i+order_max] = PARA_ARs[index_c*dim_AR1];
		for (j=1;j<dim_AR1;j++) {
			m_series[i+order_max] += (
				PARA_ARs[index_c*dim_AR1+j]*
				Y_series[i+order_max-j]);
			}
		e_series[i+order_max] = Y_series[i+order_max] - m_series[i+order_max];
		h2_series[i+order_max] = PARA_GARCHs[index_c*dim_GARCH1];
		for (j=1;j<dim_GARCH_a1;j++) {
			h2_series[i+order_max] += (
				PARA_GARCHs[index_c*dim_GARCH1+j]*
				e_series[i+order_max-j]*
				e_series[i+order_max-j]);
			}
		for (j=0;j<dim_GARCH_b1;j++) {
			h2_series[i+order_max] += (
				PARA_GARCHs[index_c*dim_GARCH1+dim_GARCH_a1+j]*
				h2_series[i+order_max-1-j]);
			}
		ln_d_norm(Y_series+i+order_max,
			m_series+i+order_max,
			h2_series+i+order_max,
			&ln_d_norm_s);
		*logLikelihood_s += ln_d_norm_s;
		}
	// for (i=0;i<*Nsize;i++) { Rprintf("%lf\t",Y_series[i+order_max]); } Rprintf("\n");
	// for (i=0;i<*Nsize;i++) { Rprintf("%lf\t",m_series[i+order_max]); } Rprintf("\n");
	// for (i=0;i<*Nsize;i++) { Rprintf("%lf\t",h2_series[i+order_max]); } Rprintf("\n");
	}

void logLikelihood(
	int *C, 				// [1]
	int *e, 				// [2]
	int *Nsize, 			// [3]
	int *np, 				// [4]
	int *orders,			// [5]  // (max.order,ar.order,garch.a.order,garch.b.order)
	int *dims,				// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARA_ARs,		// [7]  // rep(c(0,rep(0,ar.order)),n.size)
	double *PARA_GARCHs,	// [8]  // rep(c(1,rep(0,garch.a.order),rep(0,garch.b.order)),n.size)
	double *Y_series,		// [9]  // n.size + max.order
	double *m_series,		// [10] // n.size + max.order
	double *e_series,		// [11] // n.size + max.order
	double *h2_series,		// [12] // n.size + max.order
	double *logLikelihood_s	// [13]
	) {
	static int iii = 0;
	logLikelihood_i_to_n(C,e,Nsize,np,orders,dims,
		PARA_ARs,PARA_GARCHs,
		Y_series,m_series,e_series,h2_series,
		logLikelihood_s,&iii);
	}

void logcumsum(double *logA,double *Alogcumsum,int *AN) {
	int i;
	double diff;
	for (i=1;i<(*AN);i++) {
		diff = logA[i-1] - logA[i];
		if (diff>0) {
			logA[i] = logA[i-1] + log(1.0+exp(-diff));
			} else {
			logA[i] = logA[i] + log(1.0+exp(diff));
			}
		}
	*Alogcumsum = logA[*AN-1];
	for (i=0;i<*AN;i++) { logA[i] = logA[i] - *Alogcumsum; }
	}

void random_GARCH_pq(
	int *orders,
	int *dims,
	double *Random_GARCHs,
	double *PRIOR_GARCHs
	) {
	int j;
	// int dim_GARCH1 = dims[2];
	for (j=0;j<dims[2];j++) {
		r_norm(PRIOR_GARCHs+j*2,PRIOR_GARCHs+j*2+1,Random_GARCHs+j);
		if (Random_GARCHs[j]<0) { Random_GARCHs[j] *= (-1); }
		}
	}

void random_AR_p(
	int *orders,
	int *dims,
	double *Random_ARs,
	double *PRIOR_ARs
	) {
	int j;
	// int dim_AR1 = dims[0];
	for (j=0;j<dims[0];j++) {
		r_norm(PRIOR_ARs+j*2,PRIOR_ARs+j*2+1,Random_ARs+j);
		}
	}

// pred sampling step -------------------------- Begin

void logLikelihood_pred_PD(
	int *C, 				// [1]
	int *e, 				// [2]
	int *Nsize, 			// [3]
	int *np, 				// [4]
	int *orders,			// [5]  // (max.order,ar.order,garch.a.order,garch.b.order)
	int *dims,				// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARAs,			// [7]
	double *PARA_ARs,		// [8]  // rep(c(0,rep(0,ar.order)),n.size)
	double *PARA_GARCHs,	// [9]  // rep(c(1,rep(0,garch.a.order),rep(0,garch.b.order)),n.size);
	double *PRIOR_ARs,		// [10]
	double *PRIOR_GARCHs,	// [11]
	double *Y_series,		// [12] // n.size + max.order
	double *m_series,		// [13] // n.size + max.order
	double *e_series,		// [14] // n.size + max.order
	double *h2_series,		// [15] // n.size + max.order
	double *logLikelihood_s,// [16]
	double *Pred_X,			// [17]
	double *Pred_lnY,		// [18]
	int *Pred_n,			// [19]
	int *mix_ar,			// [20]
	int *mix_garch			// [21]
	) {
	int i,j,k;
	int order_max = orders[0];
	int dim_AR1 = dims[0];
	int dim_GARCH1 = dims[2];
	int dim_GARCH_a1 = dims[4];
	int dim_GARCH_b1 = dims[6];
	double Pred_m_series,Pred_h2_series;
	double ln_r_unif_s;
	double lncumsumW;
	int lengthW;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *Pred_PARA_ARs; Pred_PARA_ARs = (double *) calloc(dim_AR1,sizeof(double));
	double *Pred_PARA_GARCHs; Pred_PARA_GARCHs = (double *) calloc(dim_GARCH1,sizeof(double));

	double *PD_alpha; PD_alpha = (double *) calloc(1,sizeof(double));
	double *PD_q; PD_q = (double *) calloc(1,sizeof(double));

	*PD_alpha = *PARAs;
	*PD_q = *(PARAs+3);

	logLikelihood(C,e,Nsize,np,orders,dims,
		PARA_ARs,PARA_GARCHs,
		Y_series,m_series,e_series,h2_series,
		logLikelihood_s);

	for (j=0;j<(*np);j++) { lnW[j] = log(e[j]-(*PD_alpha)); }
	lnW[*np] = log((*PD_q)+(*PD_alpha)*(*np));
	lengthW = *np + 1;
	logcumsum(lnW,&lncumsumW,&lengthW);

	ln_r_unif(&ln_r_unif_s);
	k = *np;
	for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { k = j; j = *np; } }
	if (k<*np) {
		for (j=0;j<dim_AR1;j++) { Pred_PARA_ARs[j] = PARA_ARs[k*dim_AR1+j]; }
		for (j=0;j<dim_GARCH1;j++) { Pred_PARA_GARCHs[j] = PARA_GARCHs[k*dim_GARCH1+j]; }
		} else {
		if (*mix_ar==1) {
			random_AR_p(orders,dims,Pred_PARA_ARs,PRIOR_ARs);
			} else {
			for (j=0;j<dim_AR1;j++) { Pred_PARA_ARs[j] = PARA_ARs[j]; }
			}
		if (*mix_garch==1) {
			random_GARCH_pq(orders,dims,Pred_PARA_GARCHs,PRIOR_GARCHs);
			} else {
			for (j=0;j<dim_GARCH1;j++) { Pred_PARA_GARCHs[j] = PARA_GARCHs[j]; }
			}
		}

	Pred_m_series = Pred_PARA_ARs[0];
	for (j=1;j<dim_AR1;j++) {
		Pred_m_series += (
			Pred_PARA_ARs[j]*
			Y_series[*Nsize+order_max-j]);
		}
	Pred_h2_series = Pred_PARA_GARCHs[0];
	for (j=1;j<dim_GARCH_a1;j++) {
		Pred_h2_series += (
			Pred_PARA_GARCHs[j]*
			e_series[*Nsize+order_max-j]*
			e_series[*Nsize+order_max-j]);
		}
	for (j=0;j<dim_GARCH_b1;j++) {
		Pred_h2_series += (
			Pred_PARA_GARCHs[dim_GARCH_a1+j]*
			h2_series[*Nsize+order_max-1-j]);
		}
	for (i=0;i<*Pred_n;i++) {
		ln_d_norm(Pred_X+i,
			&Pred_m_series,
			&Pred_h2_series,
			Pred_lnY+i);
		}

	free(lnW);
	free(Pred_PARA_ARs);
	free(Pred_PARA_GARCHs);
	free(PD_alpha);
	free(PD_q);
	}

void logLikelihood_pred_NGG(
	int *C, 				// [1]
	int *e, 				// [2]
	int *Nsize, 			// [3]
	int *np, 				// [4]
	int *orders,			// [5]  // (max.order,ar.order,garch.a.order,garch.b.order)
	int *dims,				// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARAs,			// [7]
	double *PARA_ARs,		// [8]  // rep(c(0,rep(0,ar.order)),n.size)
	double *PARA_GARCHs,	// [9]  // rep(c(1,rep(0,garch.a.order),rep(0,garch.b.order)),n.size);
	double *PRIOR_ARs,		// [10]
	double *PRIOR_GARCHs,	// [11]
	double *Y_series,		// [12] // n.size + max.order
	double *m_series,		// [13] // n.size + max.order
	double *e_series,		// [14] // n.size + max.order
	double *h2_series,		// [15] // n.size + max.order
	double *logLikelihood_s,// [16]
	double *Pred_X,			// [17]
	double *Pred_lnY,		// [18]
	int *Pred_n,			// [19]
	int *mix_ar,			// [20]
	int *mix_garch			// [21]
	) {
	int i,j,k;
	int order_max = orders[0];
	int dim_AR1 = dims[0];
	int dim_GARCH1 = dims[2];
	int dim_GARCH_a1 = dims[4];
	int dim_GARCH_b1 = dims[6];
	double Pred_m_series,Pred_h2_series;
	double ln_r_unif_s;
	double gamma_a,gamma_b,r_gamma_s,logpsi,V;
	double new_lnW;
	double lncumsumW;
	int lengthW;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *Pred_PARA_ARs; Pred_PARA_ARs = (double *) calloc(dim_AR1,sizeof(double));
	double *Pred_PARA_GARCHs; Pred_PARA_GARCHs = (double *) calloc(dim_GARCH1,sizeof(double));

	double *NGG_alpha; NGG_alpha = (double *) calloc(1,sizeof(double));
	double *NGG_theta; NGG_theta = (double *) calloc(1,sizeof(double));
	double *NGG_b; NGG_b = (double *) calloc(1,sizeof(double));

	*NGG_alpha = *PARAs;
	*NGG_theta = *(PARAs+1);
	*NGG_b = *(PARAs+2);

	logLikelihood(C,e,Nsize,np,orders,dims,
		PARA_ARs,PARA_GARCHs,
		Y_series,m_series,e_series,h2_series,
		logLikelihood_s);

	// Rprintf("Reach here\n");
	// Rprintf("%10.5lf\n",*NGG_alpha);
	// Rprintf("%10.5lf\n",*NGG_theta);
	// Rprintf("%10.5lf\n",*NGG_b);

	do {
		gamma_a = *np;
		gamma_b = *NGG_theta / *NGG_alpha;
		r_gamma(&gamma_a,&gamma_b,&r_gamma_s);

		V = exp(log(r_gamma_s)/(*NGG_alpha));
		// logpsi = -((*NGG_theta)/(*NGG_alpha))
		// 	* (pow(V+*NGG_b,*NGG_alpha)-pow(V,*NGG_alpha))
		// 	+ (*Nsize-(*NGG_alpha)*(*np))*(log(V)-log(V+*NGG_b));
		logpsi = -((*NGG_theta)/(*NGG_alpha)
			*(exp((*NGG_alpha)*log(V+*NGG_b))-exp((*NGG_alpha)*log(V))))
			+ (*Nsize-(*NGG_alpha)*(*np))*(log(V)-log(V+*NGG_b));
		ln_r_unif(&ln_r_unif_s);

		// Rprintf("%10.5lf,%10.5lf,",ln_r_unif_s,logpsi);
		// Rprintf("%10.5lf,%10.5lf,%10.5lf,%10.5lf\n",V,r_gamma_s,gamma_a,gamma_b);

		} while (ln_r_unif_s>logpsi);
		
	// Rprintf("Reach here\n");
		
		
	new_lnW = log(*NGG_theta) + (*NGG_alpha) * log(V+(*NGG_b));

	for (j=0;j<(*np);j++) {
		lnW[j] = log(e[j]-(*NGG_alpha))
				- log((*Nsize)-(*np)*(*NGG_alpha)+exp(new_lnW));
		}
	lnW[*np] = new_lnW - log((*Nsize)-(*np)*(*NGG_alpha)+exp(new_lnW));
	lengthW = *np + 1;
	logcumsum(lnW,&lncumsumW,&lengthW);

	ln_r_unif(&ln_r_unif_s);
	k = *np;
	for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { k = j; j = *np; } }
	if (k<*np) {
		for (j=0;j<dim_AR1;j++) { Pred_PARA_ARs[j] = PARA_ARs[k*dim_AR1+j]; }
		for (j=0;j<dim_GARCH1;j++) { Pred_PARA_GARCHs[j] = PARA_GARCHs[k*dim_GARCH1+j]; }
		} else {
		if (*mix_ar==1) {
			random_AR_p(orders,dims,Pred_PARA_ARs,PRIOR_ARs);
			} else {
			for (j=0;j<dim_AR1;j++) { Pred_PARA_ARs[j] = PARA_ARs[j]; }
			}
		if (*mix_garch==1) {
			random_GARCH_pq(orders,dims,Pred_PARA_GARCHs,PRIOR_GARCHs);
			} else {
			for (j=0;j<dim_GARCH1;j++) { Pred_PARA_GARCHs[j] = PARA_GARCHs[j]; }
			}
		}
		
	Pred_m_series = Pred_PARA_ARs[0];
	for (j=1;j<dim_AR1;j++) {
		Pred_m_series += (
			Pred_PARA_ARs[j]*
			Y_series[*Nsize+order_max-j]);
		}
	Pred_h2_series = Pred_PARA_GARCHs[0];
	for (j=1;j<dim_GARCH_a1;j++) {
		Pred_h2_series += (
			Pred_PARA_GARCHs[j]*
			e_series[*Nsize+order_max-j]*
			e_series[*Nsize+order_max-j]);
		}
	for (j=0;j<dim_GARCH_b1;j++) {
		Pred_h2_series += (
			Pred_PARA_GARCHs[dim_GARCH_a1+j]*
			h2_series[*Nsize+order_max-1-j]);
		}
	for (i=0;i<*Pred_n;i++) {
		ln_d_norm(Pred_X+i,
			&Pred_m_series,
			&Pred_h2_series,
			Pred_lnY+i);
		}

	free(lnW);
	free(Pred_PARA_ARs);
	free(Pred_PARA_GARCHs);
	free(NGG_alpha);
	free(NGG_theta);
	free(NGG_b);
	}

void logLikelihood_pred_NON(
	int *C, 				// [1]
	int *e, 				// [2]
	int *Nsize, 			// [3]
	int *np, 				// [4]
	int *orders,			// [5]  // (max.order,ar.order,garch.a.order,garch.b.order)
	int *dims,				// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARAs,			// [7]
	double *PARA_ARs,		// [8]  // rep(c(0,rep(0,ar.order)),n.size)
	double *PARA_GARCHs,	// [9]  // rep(c(1,rep(0,garch.a.order),rep(0,garch.b.order)),n.size);
	double *PRIOR_ARs,		// [10]
	double *PRIOR_GARCHs,	// [11]
	double *Y_series,		// [12] // n.size + max.order
	double *m_series,		// [13] // n.size + max.order
	double *e_series,		// [14] // n.size + max.order
	double *h2_series,		// [15] // n.size + max.order
	double *logLikelihood_s,// [16]
	double *Pred_X,			// [17]
	double *Pred_lnY,		// [18]
	int *Pred_n,			// [19]
	int *mix_ar,			// [20]
	int *mix_garch			// [21]
	) {
	int i,j;
	int order_max = orders[0];
	int dim_AR1 = dims[0];
	int dim_GARCH1 = dims[2];
	int dim_GARCH_a1 = dims[4];
	int dim_GARCH_b1 = dims[6];
	double Pred_m_series,Pred_h2_series;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *Pred_PARA_ARs; Pred_PARA_ARs = (double *) calloc(dim_AR1,sizeof(double));
	double *Pred_PARA_GARCHs; Pred_PARA_GARCHs = (double *) calloc(dim_GARCH1,sizeof(double));

	double *NGG_alpha; NGG_alpha = (double *) calloc(1,sizeof(double));
	double *NGG_theta; NGG_theta = (double *) calloc(1,sizeof(double));
	double *NGG_b; NGG_b = (double *) calloc(1,sizeof(double));

	*NGG_alpha = *PARAs;
	*NGG_theta = *(PARAs+1);
	*NGG_b = *(PARAs+2);

	logLikelihood(C,e,Nsize,np,orders,dims,
		PARA_ARs,PARA_GARCHs,
		Y_series,m_series,e_series,h2_series,
		logLikelihood_s);

	// Rprintf("Reach here\n");
	// Rprintf("%10.5lf\n",*NGG_alpha);
	// Rprintf("%10.5lf\n",*NGG_theta);
	// Rprintf("%10.5lf\n",*NGG_b);

	for (j=0;j<dim_AR1;j++) { Pred_PARA_ARs[j] = PARA_ARs[j]; }
	for (j=0;j<dim_GARCH1;j++) { Pred_PARA_GARCHs[j] = PARA_GARCHs[j]; }
		
	Pred_m_series = Pred_PARA_ARs[0];
	for (j=1;j<dim_AR1;j++) {
		Pred_m_series += (
			Pred_PARA_ARs[j]*
			Y_series[*Nsize+order_max-j]);
		}
	Pred_h2_series = Pred_PARA_GARCHs[0];
	for (j=1;j<dim_GARCH_a1;j++) {
		Pred_h2_series += (
			Pred_PARA_GARCHs[j]*
			e_series[*Nsize+order_max-j]*
			e_series[*Nsize+order_max-j]);
		}
	for (j=0;j<dim_GARCH_b1;j++) {
		Pred_h2_series += (
			Pred_PARA_GARCHs[dim_GARCH_a1+j]*
			h2_series[*Nsize+order_max-1-j]);
		}
	for (i=0;i<*Pred_n;i++) {
		ln_d_norm(Pred_X+i,
			&Pred_m_series,
			&Pred_h2_series,
			Pred_lnY+i);
		}

	free(lnW);
	free(Pred_PARA_ARs);
	free(Pred_PARA_GARCHs);
	free(NGG_alpha);
	free(NGG_theta);
	free(NGG_b);
	}

// pred sampling step -------------------------- End

// partition sampling step --------------------- Begin

void sample_PD_GIB_1STEP_MGARCH(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARAs,					// [7]
	double *CURVALUE_PARA_ARs,		// [8]
	double *PROPOSAL_PARA_ARs,		// [9]
	double *CURVALUE_PARA_GARCHs,	// [10]
	double *PROPOSAL_PARA_GARCHs,	// [11]
	double *PRIOR_ARs,				// [12]
	double *PRIOR_GARCHs,			// [13]
	double *Y_series,				// [14]
	double *CURVALUE_m_series,		// [15]
	double *PROPOSAL_m_series,		// [16]
	double *CURVALUE_e_series,		// [17]
	double *PROPOSAL_e_series,		// [18]
	double *CURVALUE_h2_series,		// [19]
	double *PROPOSAL_h2_series,		// [20]
	int *Option,					// [21]
	int *max_R,						// [22]
	int *approx_N,					// [23]
	int *mix_ar,					// [24]
	int *mix_garch					// [25]
	) {
	int i,j,k,index_c;
	double ln_r_unif_s;
	double approx_lnW;
	// double diff_lnW;
	double PROPOSAL_logLikelihood_s;
	int dim_AR1 = dims[0];
	int dim_GARCH1 = dims[2];
	int start_i;
	double lncumsumW;
	int lengthW;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *PD_alpha; PD_alpha = (double *) calloc(1,sizeof(double));
	double *PD_q; PD_q = (double *) calloc(1,sizeof(double));

	double *Random_PARA_ARs; Random_PARA_ARs = (double *) calloc((*approx_N)*dim_AR1,sizeof(double));
	double *Random_PARA_GARCHs; Random_PARA_GARCHs = (double *) calloc((*approx_N)*dim_GARCH1,sizeof(double));
	double *Random_lnW; Random_lnW = (double *) calloc(*approx_N,sizeof(double));

	*PD_alpha = *PARAs;
	*PD_q = *(PARAs+3);

	if (*mix_ar==1||*mix_garch==1) {

		for (i=0;i<*Nsize;i++) {

			index_c = C[i] - 1;
			e[index_c] = e[index_c] - 1;
			if (e[index_c]<1) {
				for (j=0;j<(*Nsize);j++) {
					if (C[j]==(*np)) { C[j] = C[i]; }
					}
				*np = *np - 1;
				e[index_c] = e[*np];
				for (j=0;j<dim_GARCH1;j++) {
					CURVALUE_PARA_GARCHs[index_c*dim_GARCH1+j] =
					CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+j];
					}
				for (j=0;j<dim_AR1;j++) {
					CURVALUE_PARA_ARs[index_c*dim_AR1+j] =
					CURVALUE_PARA_ARs[(*np)*dim_AR1+j];
					}
				}
			e[*np] = 0;

			logLikelihood(C,e,Nsize,np,orders,dims,
				CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);

			start_i = i;
			for (j=0;j<(*np);j++) {
				C[i] = j + 1;
				logLikelihood_i_to_n(C,e,Nsize,np,orders,dims,
					CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
					Y_series,
					PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
					&PROPOSAL_logLikelihood_s,&start_i);
				lnW[j] = log(e[j]-(*PD_alpha)) + PROPOSAL_logLikelihood_s;
				}

			C[i] = *np + 1;
			approx_lnW = 0.0;
			// diff_lnW = 0.0;
			for (j=0;j<*approx_N;j++) {
				if (*mix_ar==1) {
					random_AR_p(orders,dims,Random_PARA_ARs+j*dim_AR1,PRIOR_ARs);
					for (k=0;k<dim_AR1;k++) {
						CURVALUE_PARA_ARs[(*np)*dim_AR1+k] =
						Random_PARA_ARs[j*dim_AR1+k];
						}
					} else {
					for (k=0;k<dim_AR1;k++) {
						CURVALUE_PARA_ARs[(*np)*dim_AR1+k] =
						CURVALUE_PARA_ARs[k];
						}
					}
				if (*mix_garch==1) {
					random_GARCH_pq(orders,dims,Random_PARA_GARCHs+j*dim_GARCH1,PRIOR_GARCHs);
					for (k=0;k<dim_GARCH1;k++) {
						CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+k] =
						Random_PARA_GARCHs[j*dim_GARCH1+k];
						}
					} else {
					for (k=0;k<dim_GARCH1;k++) {
						CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+k] =
						CURVALUE_PARA_GARCHs[k];
						}
					}
				logLikelihood_i_to_n(C,e,Nsize,np,orders,dims,
					CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
					Y_series,
					PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
					&PROPOSAL_logLikelihood_s,&start_i);
				Random_lnW[j] = PROPOSAL_logLikelihood_s;
				}
			logcumsum(Random_lnW,&approx_lnW,approx_N);

			lnW[*np] = log((*PD_q)+(*PD_alpha)*(*np)) + approx_lnW - log(*approx_N);
			lengthW = *np + 1;
			logcumsum(lnW,&lncumsumW,&lengthW);

			ln_r_unif(&ln_r_unif_s);
			index_c = *np;
			for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { index_c = j; j = *np; } }
			C[i] = index_c + 1;
			e[index_c] = e[index_c] + 1;

			if (index_c==(*np)) {
				ln_r_unif(&ln_r_unif_s);
				k = *approx_N - 1;
				for (j=0;j<(*approx_N-1);j++) {
					if (ln_r_unif_s<Random_lnW[j]) { k = j; j = *approx_N; }
					}
				if (*mix_ar==1) {
					for (j=0;j<dim_AR1;j++) {
						CURVALUE_PARA_ARs[(*np)*dim_AR1+j] =
						Random_PARA_ARs[k*dim_AR1+j];
						}
					}
				if (*mix_garch==1) {
					for (j=0;j<dim_GARCH1;j++) {
						CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+j] =
						Random_PARA_GARCHs[k*dim_GARCH1+j];
						}
					}
				*np = *np + 1;
				}

			// for (j=0;j<*approx_N;j++) { Rprintf("%lf\t",Random_lnW[j]); } Rprintf("\n");
			// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",C[j]); } Rprintf("\n");
			// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",e[j]); } Rprintf("\n");
			// for (j=0;j<*Nsize;j++) { Rprintf("%lf\t",lnW[j]); } Rprintf("\n");

			}
		}

	free(Random_PARA_ARs);
	free(Random_PARA_GARCHs);
	free(Random_lnW);
	free(PD_alpha);
	free(PD_q);
	free(lnW);
	}

void sample_NON_GIB_1STEP_MGARCH(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARAs,					// [7]
	double *CURVALUE_PARA_ARs,		// [8]
	double *PROPOSAL_PARA_ARs,		// [9]
	double *CURVALUE_PARA_GARCHs,	// [10]
	double *PROPOSAL_PARA_GARCHs,	// [11]
	double *PRIOR_ARs,				// [12]
	double *PRIOR_GARCHs,			// [13]
	double *Y_series,				// [14]
	double *CURVALUE_m_series,		// [15]
	double *PROPOSAL_m_series,		// [16]
	double *CURVALUE_e_series,		// [17]
	double *PROPOSAL_e_series,		// [18]
	double *CURVALUE_h2_series,		// [19]
	double *PROPOSAL_h2_series,		// [20]
	int *Option,					// [21]
	int *max_R,						// [22]
	int *approx_N,					// [23]
	int *mix_ar,					// [24]
	int *mix_garch					// [25]
	) {

	}

void sample_NGG_GIB_2STEP_MGARCH(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *PARAs,					// [7]
	double *CURVALUE_PARA_ARs,		// [8]
	double *PROPOSAL_PARA_ARs,		// [9]
	double *CURVALUE_PARA_GARCHs,	// [10]
	double *PROPOSAL_PARA_GARCHs,	// [11]
	double *PRIOR_ARs,				// [12]
	double *PRIOR_GARCHs,			// [13]
	double *Y_series,				// [14]
	double *CURVALUE_m_series,		// [15]
	double *PROPOSAL_m_series,		// [16]
	double *CURVALUE_e_series,		// [17]
	double *PROPOSAL_e_series,		// [18]
	double *CURVALUE_h2_series,		// [19]
	double *PROPOSAL_h2_series,		// [20]
	int *Option,					// [21]
	int *max_R,						// [22]
	int *approx_N,					// [23]
	int *mix_ar,					// [24]
	int *mix_garch					// [25]
	) {
	int i,j,k,index_c;
	double ln_r_unif_s;
	double approx_lnW;
	// double diff_lnW;
	double PROPOSAL_logLikelihood_s;
	double gamma_a,gamma_b,r_gamma_s,logpsi,V;
	int dim_AR1 = dims[0];
	int dim_GARCH1 = dims[2];
	int start_i;
	double lncumsumW;
	int lengthW;

	double *lnW; lnW = (double *) calloc(*Nsize,sizeof(double));
	double *NGG_alpha; NGG_alpha = (double *) calloc(1,sizeof(double));
	double *NGG_theta; NGG_theta = (double *) calloc(1,sizeof(double));
	double *NGG_b; NGG_b = (double *) calloc(1,sizeof(double));

	double *Random_PARA_ARs; Random_PARA_ARs = (double *) calloc((*approx_N)*dim_AR1,sizeof(double));
	double *Random_PARA_GARCHs; Random_PARA_GARCHs = (double *) calloc((*approx_N)*dim_GARCH1,sizeof(double));
	double *Random_lnW; Random_lnW = (double *) calloc(*approx_N,sizeof(double));

	*NGG_alpha = *PARAs;
	*NGG_theta = *(PARAs+1);
	*NGG_b = *(PARAs+2);

	if (*mix_ar==1||*mix_garch==1) {

		for (i=0;i<*Nsize;i++) {

			do {
				gamma_a = *np;
				gamma_b = *NGG_theta / *NGG_alpha;
				r_gamma(&gamma_a,&gamma_b,&r_gamma_s);
				V = exp(log(r_gamma_s)/(*NGG_alpha));
				logpsi = -((*NGG_theta)/(*NGG_alpha)) 
					* (pow(V+*NGG_b,*NGG_alpha)-pow(V,*NGG_alpha))
					+ (*Nsize-(*NGG_alpha)*(*np))*(log(V)-log(V+*NGG_b));
				ln_r_unif(&ln_r_unif_s);
				} while (ln_r_unif_s>logpsi); // U = V;

			index_c = C[i] - 1;
			e[index_c] = e[index_c] - 1;
			if (e[index_c]<1) {
				for (j=0;j<(*Nsize);j++) { if (C[j]==(*np)) { C[j] = C[i]; } }
				*np = *np - 1;
				e[index_c] = e[*np];
				for (j=0;j<dim_GARCH1;j++) {
					CURVALUE_PARA_GARCHs[index_c*dim_GARCH1+j] =
					CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+j];
					}
				for (j=0;j<dim_AR1;j++) {
					CURVALUE_PARA_ARs[index_c*dim_AR1+j] =
					CURVALUE_PARA_ARs[(*np)*dim_AR1+j];
					}
				}
			e[*np] = 0;

			logLikelihood(C,e,Nsize,np,orders,dims,
				CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);
				
			start_i = i;
			for (j=0;j<(*np);j++) {
				C[i] = j + 1;
				logLikelihood_i_to_n(C,e,Nsize,np,orders,dims,
					CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
					Y_series,
					PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
					&PROPOSAL_logLikelihood_s,&start_i);
				lnW[j] = log(e[j]-(*NGG_alpha)) + PROPOSAL_logLikelihood_s;
				}

			C[i] = *np + 1;
			approx_lnW = 0.0;
			// diff_lnW = 0.0;
			for (j=0;j<*approx_N;j++) {
				if (*mix_ar==1) {
					random_AR_p(orders,dims,Random_PARA_ARs+j*dim_AR1,PRIOR_ARs);
					for (k=0;k<dim_AR1;k++) {
						CURVALUE_PARA_ARs[(*np)*dim_AR1+k] =
						Random_PARA_ARs[j*dim_AR1+k];
						}
					} else {
					for (k=0;k<dim_AR1;k++) {
						CURVALUE_PARA_ARs[(*np)*dim_AR1+k] =
						CURVALUE_PARA_ARs[k];
						}
					}
				if (*mix_garch==1) {
					random_GARCH_pq(orders,dims,Random_PARA_GARCHs+j*dim_GARCH1,PRIOR_GARCHs);
					for (k=0;k<dim_GARCH1;k++) {
						CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+k] =
						Random_PARA_GARCHs[j*dim_GARCH1+k];
						}
					} else {
					for (k=0;k<dim_GARCH1;k++) {
						CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+k] =
						CURVALUE_PARA_GARCHs[k];
						}
					}
				logLikelihood_i_to_n(C,e,Nsize,np,orders,dims,
					CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
					Y_series,
					PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
					&PROPOSAL_logLikelihood_s,&start_i);
				Random_lnW[j] = PROPOSAL_logLikelihood_s;
				}
			logcumsum(Random_lnW,&approx_lnW,approx_N);

			lnW[*np] = log(*NGG_theta) + (*NGG_alpha) * log(V+(*NGG_b))
						+ approx_lnW - log(*approx_N);
			lengthW = *np + 1;
			logcumsum(lnW,&lncumsumW,&lengthW);

			ln_r_unif(&ln_r_unif_s);
			index_c = *np;
			for (j=0;j<(*np);j++) { if (ln_r_unif_s<lnW[j]) { index_c = j; j = *np; } }
			C[i] = index_c + 1;
			e[index_c] = e[index_c] + 1;

			if (index_c==(*np)) {
				ln_r_unif(&ln_r_unif_s);
				k = *approx_N - 1;
				for (j=0;j<(*approx_N-1);j++) {
					if (ln_r_unif_s<Random_lnW[j]) { k = j; j = *approx_N; }
					}
				if (*mix_ar==1) {
					for (j=0;j<dim_AR1;j++) {
						CURVALUE_PARA_ARs[(*np)*dim_AR1+j] =
						Random_PARA_ARs[k*dim_AR1+j];
						}
					}
				if (*mix_garch==1) {
					for (j=0;j<dim_GARCH1;j++) {
						CURVALUE_PARA_GARCHs[(*np)*dim_GARCH1+j] =
						Random_PARA_GARCHs[k*dim_GARCH1+j];
						}
					}
				*np = *np + 1;
				}

			// for (j=0;j<*approx_N;j++) { Rprintf("%lf\t",Random_lnW[j]); } Rprintf("\n");
			// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",C[j]); } Rprintf("\n");
			// for (j=0;j<*Nsize;j++) { Rprintf("%d\t",e[j]); } Rprintf("\n");
			// for (j=0;j<*Nsize;j++) { Rprintf("%lf\t",lnW[j]); } Rprintf("\n");
			}
		}

	free(Random_PARA_ARs);
	free(Random_PARA_GARCHs);
	free(Random_lnW);
	free(NGG_alpha);
	free(NGG_theta);
	free(NGG_b);
	free(lnW);
	}

// partition sampling step --------------------- End

// garch parameters sampling step -------------- Begin
	
void sample_GARCH_AMDR_MCMC(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *CURVALUE_PARA_ARs,		// [7]
	double *PROPOSAL_PARA_ARs,		// [8]
	double *CURVALUE_PARA_GARCHs,	// [9]
	double *PROPOSAL_PARA_GARCHs,	// [10]
	double *PRIOR_ARs,				// [11]
	double *PRIOR_GARCHs,			// [12]
	double *Y_series,				// [13]
	double *CURVALUE_m_series,		// [14]
	double *PROPOSAL_m_series,		// [15]
	double *CURVALUE_e_series,		// [16]
	double *PROPOSAL_e_series,		// [17]
	double *CURVALUE_h2_series,		// [18]
	double *PROPOSAL_h2_series,		// [19]
	int *Option,					// [20]
	double *AMDR_AVE_ARs,			// [21]
	double *AMDR_COV_ARs,			// [22]
	double *AMDR_AVE_GARCHs,		// [23]
	double *AMDR_COV_GARCHs,		// [24]
	int *Iter,						// [25]
	int *max_R						// [26]
	) {
	int i,j,k,l;
	double CURVALUE_logLikelihood_s,CURVALUE_logPrior_s,CURVALUE_logProposal_s;
	double PROPOSAL_logLikelihood_s,PROPOSAL_logPrior_s,PROPOSAL_logProposal_s;
	double kappa_CURVALUE,kappa_PROPOSAL,kappa_MAXIMUM;
	double acc_prob;
	double r_unif_s;
	int dim_GARCH1 = dims[2];
	int dim_GARCH2 = dims[3];
	int first_cal,index,positive,count;
	
	double *Local_AVE; Local_AVE = (double *) calloc(dim_GARCH1,sizeof(double));
	double *Local_COV; Local_COV = (double *) calloc(dim_GARCH2,sizeof(double));
	double *PRIOR_AVE; PRIOR_AVE = (double *) calloc(dim_GARCH1,sizeof(double));
	double *PRIOR_COV; PRIOR_COV = (double *) calloc(dim_GARCH2,sizeof(double));
	
	for (i=0;i<(*np);i++) {
		for (j=0;j<dim_GARCH1;j++) {
			PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
			CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
			}
		}

	for (k=0;k<dim_GARCH1;k++) { PRIOR_AVE[k] = PRIOR_GARCHs[2*k]; }
	for (k=0;k<dim_GARCH2;k++) { PRIOR_COV[k] = 0.0; }
	for (k=0;k<dim_GARCH1;k++) { PRIOR_COV[k*dim_GARCH1+k] = PRIOR_GARCHs[2*k+1]; }

	for (i=0;i<(*np);i++) {

		for (j=0;j<dim_GARCH1;j++) {
			PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
			CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
			}

		// for (k=0;k<dim_GARCH1;k++) { Local_AVE[k] = CURVALUE_PARA_GARCHs[i*dim_GARCH1+k]; }
		for (k=0;k<dim_GARCH1;k++) { Local_AVE[k] = AMDR_AVE_GARCHs[i*dim_GARCH1+k]; }
		for (k=0;k<dim_GARCH2;k++) { Local_COV[k] = 2.34 * AMDR_COV_GARCHs[i*dim_GARCH2+k]; }

		logLikelihood(C,e,Nsize,np,orders,dims,
			CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
			Y_series,
			CURVALUE_m_series,CURVALUE_e_series,CURVALUE_h2_series,
			&CURVALUE_logLikelihood_s);
		ln_d_mv_norm(CURVALUE_PARA_GARCHs+i*dim_GARCH1,Local_AVE,Local_COV,&dim_GARCH1,&CURVALUE_logProposal_s);
		ln_d_mv_norm(CURVALUE_PARA_GARCHs+i*dim_GARCH1,PRIOR_AVE,PRIOR_COV,&dim_GARCH1,&CURVALUE_logPrior_s);
		kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s - CURVALUE_logProposal_s;
		// kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s;

		kappa_MAXIMUM = 0.0;
		first_cal = 1;
		index = 0;
		positive = 0;
		count = 0;

		do {

			do {

				r_mv_norm(Local_AVE,Local_COV,&dim_GARCH1,PROPOSAL_PARA_GARCHs+i*dim_GARCH1);

				positive = 1;
				for (j=0;j<dim_GARCH1;j++) {
					if (PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j]<0.0) { positive = 0; }
					}

				} while(positive==0);

			logLikelihood(C,e,Nsize,np,orders,dims,
				CURVALUE_PARA_ARs,PROPOSAL_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);
			ln_d_mv_norm(PROPOSAL_PARA_GARCHs+i*dim_GARCH1,Local_AVE,Local_COV,&dim_GARCH1,&PROPOSAL_logProposal_s);
			ln_d_mv_norm(PROPOSAL_PARA_GARCHs+i*dim_GARCH1,PRIOR_AVE,PRIOR_COV,&dim_GARCH1,&PROPOSAL_logPrior_s);
			kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s - PROPOSAL_logProposal_s;
			// kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s;

			if (first_cal==1) {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else {
					acc_prob = exp(kappa_PROPOSAL-kappa_CURVALUE);
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					first_cal = 0;
					}
				} else {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else if (kappa_PROPOSAL<kappa_MAXIMUM) {
					acc_prob = 0.0;
					index = 0;
					} else {
					acc_prob = 1.0
						- (1.0-exp(kappa_PROPOSAL-kappa_CURVALUE))
						/ (1.0-exp(kappa_MAXIMUM-kappa_CURVALUE));
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					}
				}
			
			if (index==1) {
				for (j=0;j<dim_GARCH1;j++) {
					CURVALUE_PARA_GARCHs[i*dim_GARCH1+j] =
					PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j];
					}
				} else {
				for (j=0;j<dim_GARCH1;j++) {
					PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
					CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
					}
				}

			count++; if (count>*max_R) { index = 1; }
			} while(index==0);

		for (k=0;k<dim_GARCH1;k++) {
			AMDR_AVE_GARCHs[i*dim_GARCH1+k]
				= AMDR_AVE_GARCHs[i*dim_GARCH1+k]
				+ (CURVALUE_PARA_GARCHs[i*dim_GARCH1+k]-AMDR_AVE_GARCHs[i*dim_GARCH1+k])
				/ (*Iter+1);
			}

		for (k=0;k<dim_GARCH1;k++) {
			for (l=0;l<dim_GARCH1;l++) {
				AMDR_COV_GARCHs[i*dim_GARCH2+k*dim_GARCH1+l]
					= AMDR_COV_GARCHs[i*dim_GARCH2+k*dim_GARCH1+l] +
					((CURVALUE_PARA_GARCHs[i*dim_GARCH1+k]-AMDR_AVE_GARCHs[i*dim_GARCH1+k])
					*(CURVALUE_PARA_GARCHs[i*dim_GARCH1+l]-AMDR_AVE_GARCHs[i*dim_GARCH1+l])
					-AMDR_COV_GARCHs[i*dim_GARCH2+k*dim_GARCH1+l]) / (*Iter+1);
				}
			}
		}

	free(Local_AVE);
	free(Local_COV);
	free(PRIOR_AVE);
	free(PRIOR_COV);
	}

void sample_GARCH_pq_MCMC(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *CURVALUE_PARA_ARs,		// [7]
	double *PROPOSAL_PARA_ARs,		// [8]
	double *CURVALUE_PARA_GARCHs,	// [9]
	double *PROPOSAL_PARA_GARCHs,	// [10]
	double *PRIOR_ARs,				// [11]
	double *PRIOR_GARCHs,			// [12]
	double *Y_series,				// [13]
	double *CURVALUE_m_series,		// [14]
	double *PROPOSAL_m_series,		// [15]
	double *CURVALUE_e_series,		// [16]
	double *PROPOSAL_e_series,		// [17]
	double *CURVALUE_h2_series,		// [18]
	double *PROPOSAL_h2_series,		// [19]
	int *Option,					// [20]
	double *AMDR_AVE_ARs,			// [21]
	double *AMDR_COV_ARs,			// [22]
	double *AMDR_AVE_GARCHs,		// [23]
	double *AMDR_COV_GARCHs,		// [24]
	int *Iter,						// [25]
	int *max_R						// [26]
	) {
	int i,j,k,l;
	int index_c;
	double CURVALUE_logLikelihood_s,CURVALUE_logPrior_s,CURVALUE_logProposal_s;
	double PROPOSAL_logLikelihood_s,PROPOSAL_logPrior_s,PROPOSAL_logProposal_s;
	double kappa_CURVALUE,kappa_PROPOSAL,kappa_MAXIMUM;
	double acc_prob;
	double r_unif_s;
	int order_max = orders[0];
	int dim_GARCH1 = dims[2];
	int dim_GARCH2 = dims[3];
	int dim_GARCH_a1 = dims[4];
	int dim_GARCH_a2 = dims[5];
	int dim_GARCH_b1 = dims[6];
	int dim_GARCH_b2 = dims[7];
	int first_cal,index,positive,count;

	double *e2_tilde; e2_tilde = calloc((*Nsize)+order_max,sizeof(double));
	double *i2_tilde; i2_tilde = calloc((*Nsize)+order_max,sizeof(double));
	double *e2; e2 = calloc((*Nsize)+order_max,sizeof(double));
	double *h4; h4 = calloc((*Nsize)+order_max,sizeof(double));

	double *XY; XY = calloc(dim_GARCH1,sizeof(double));
	double *XX; XX = calloc(dim_GARCH2,sizeof(double));
	double *Yi; Yi = calloc(1,sizeof(double));
	double *Xi; Xi = calloc(dim_GARCH1,sizeof(double));
	double *Zi; Zi = calloc(((*Nsize)+order_max)*dim_GARCH_b1,sizeof(double));
	double *inv_XX; inv_XX = calloc(dim_GARCH2,sizeof(double));
	double *Local_AVE; Local_AVE = (double *) calloc(dim_GARCH1,sizeof(double));
	double *Local_COV; Local_COV = (double *) calloc(dim_GARCH2,sizeof(double));
	double *PRIOR_AVE; PRIOR_AVE = (double *) calloc(dim_GARCH1,sizeof(double));
	double *PRIOR_COV; PRIOR_COV = (double *) calloc(dim_GARCH2,sizeof(double));
	double *PRIOR_inv_COV; PRIOR_inv_COV = (double *) calloc(dim_GARCH2,sizeof(double));

	for (i=0;i<(*np);i++) {
		for (j=0;j<dim_GARCH1;j++) {
			PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
			CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
			}
		}

	for (k=0;k<dim_GARCH_a1;k++) { PRIOR_AVE[k] = PRIOR_GARCHs[2*k]; }
	for (k=0;k<dim_GARCH_a2;k++) { PRIOR_COV[k] = 0.0; }
	for (k=0;k<dim_GARCH_a1;k++) { PRIOR_COV[k*dim_GARCH_a1+k] = PRIOR_GARCHs[2*k+1]; }

	for (i=0;i<(*np);i++) {

		for (j=0;j<dim_GARCH1;j++) {
			PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] = 
			CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
			}

		logLikelihood(C,e,Nsize,np,orders,dims,
			CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
			Y_series,
			CURVALUE_m_series,CURVALUE_e_series,CURVALUE_h2_series,
			&CURVALUE_logLikelihood_s);

		inv(PRIOR_COV,&dim_GARCH_a1,PRIOR_inv_COV);
		for (k=0;k<dim_GARCH_a1;k++) {
			XY[k] = 0.0;
			for (l=0;l<dim_GARCH_a1;l++) {
				XY[k] += (PRIOR_inv_COV[k*dim_GARCH_a1+l]*PRIOR_AVE[l]);
				XX[k*dim_GARCH_a1+l] = PRIOR_inv_COV[k*dim_GARCH_a1+l];
				}
			}

		for (j=0;j<(*Nsize);j++) {
			e2[j+order_max] = CURVALUE_e_series[j+order_max] * CURVALUE_e_series[j+order_max];
			h4[j+order_max] = CURVALUE_h2_series[j+order_max] * CURVALUE_h2_series[j+order_max];
			index_c = C[j] - 1;
			e2_tilde[j+order_max] = e2[j+order_max];
			i2_tilde[j+order_max] = 1.0;
			for (k=0;k<dim_GARCH_b1;k++) {
				e2_tilde[j+order_max] += (CURVALUE_PARA_GARCHs[index_c*dim_GARCH1+dim_GARCH_a1+k]*e2_tilde[j+order_max-1-k]);
				i2_tilde[j+order_max] += (CURVALUE_PARA_GARCHs[index_c*dim_GARCH1+dim_GARCH_a1+k]*i2_tilde[j+order_max-1-k]);
				}
			if (index_c==i) {
				Yi[0] = e2[j+order_max];
				Xi[0] = i2_tilde[j+order_max];
				for (k=1;k<dim_GARCH_a1;k++) {
					Xi[k] = e2_tilde[j+order_max-k];
					}
				for (k=0;k<dim_GARCH_a1;k++) {
					XY[k] += (Xi[k]*Yi[0]*0.5/h4[j+order_max]);
					for (l=0;l<dim_GARCH_a1;l++) {
						XX[k*dim_GARCH_a1+l] += (Xi[k]*Xi[l]*0.5/h4[j+order_max]);
						}
					}
				}
			}

		inv(XX,&dim_GARCH_a1,inv_XX);
		for (k=0;k<dim_GARCH_a1;k++) {
			Local_AVE[k] = 0.0;
			for (l=0;l<dim_GARCH_a1;l++) {
				Local_AVE[k] += (inv_XX[k*dim_GARCH_a1+l]*XY[l]);
				Local_COV[k*dim_GARCH_a1+l] = 10.0 * inv_XX[k*dim_GARCH_a1+l];
				// if (k==l) Local_COV[k*dim_GARCH_a1+l] = 0.1; else Local_COV[k*dim_GARCH_a1+l] = 0.0;
				}
			}

		ln_d_mv_norm(CURVALUE_PARA_GARCHs+i*dim_GARCH1,Local_AVE,Local_COV,&dim_GARCH_a1,&CURVALUE_logProposal_s);
		ln_d_mv_norm(CURVALUE_PARA_GARCHs+i*dim_GARCH1,PRIOR_AVE,PRIOR_COV,&dim_GARCH_a1,&CURVALUE_logPrior_s);
		kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s - CURVALUE_logProposal_s;

		kappa_MAXIMUM = 0.0;
		first_cal = 1;
		index = 0;
		positive = 0;
		count = 0;

		do {

			do {
				r_mv_norm(Local_AVE,Local_COV,&dim_GARCH_a1,PROPOSAL_PARA_GARCHs+i*dim_GARCH1);
				positive = 1;
				for (j=0;j<dim_GARCH_a1;j++) {
					if (PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j]<0.0) { positive = 0; }
					}
				} while(positive==0);

			logLikelihood(C,e,Nsize,np,orders,dims,
				CURVALUE_PARA_ARs,PROPOSAL_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);
			ln_d_mv_norm(PROPOSAL_PARA_GARCHs+i*dim_GARCH1,Local_AVE,Local_COV,&dim_GARCH_a1,&PROPOSAL_logProposal_s);
			ln_d_mv_norm(PROPOSAL_PARA_GARCHs+i*dim_GARCH1,PRIOR_AVE,PRIOR_COV,&dim_GARCH_a1,&PROPOSAL_logPrior_s);
			kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s - PROPOSAL_logProposal_s;

			if (first_cal==1) {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else {
					acc_prob = exp(kappa_PROPOSAL-kappa_CURVALUE);
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					first_cal = 0;
					}
				} else {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else if (kappa_PROPOSAL<kappa_MAXIMUM) {
					acc_prob = 0.0;
					index = 0;
					} else {
					acc_prob = 1.0
						- (1.0-exp(kappa_PROPOSAL-kappa_CURVALUE))
						/ (1.0-exp(kappa_MAXIMUM-kappa_CURVALUE));
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					}
				}
			
			if (index==1) {
				for (j=0;j<dim_GARCH_a1;j++) {
					CURVALUE_PARA_GARCHs[i*dim_GARCH1+j] =
					PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j];
					}
				} else {
				for (j=0;j<dim_GARCH_a1;j++) {
					PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
					CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
					}
				}

			count++; if (count>*max_R) { index = 1; }
			} while(index==0);

		}

	for (k=0;k<dim_GARCH_b1;k++) { PRIOR_AVE[k] = PRIOR_GARCHs[2*(dim_GARCH_a1+k)]; }
	for (k=0;k<dim_GARCH_b2;k++) { PRIOR_COV[k] = 0.0; }
	for (k=0;k<dim_GARCH_b1;k++) { PRIOR_COV[k*dim_GARCH_b1+k] = PRIOR_GARCHs[2*(dim_GARCH_a1+k)+1]; }

	for (i=0;i<(*np);i++) {

		for (j=0;j<dim_GARCH1;j++) {
			PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
			CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
			}

		// for (j=0;j<dim_GARCH1;j++) { Rprintf("%lf\t",CURVALUE_PARA_GARCHs[i*dim_GARCH1+j]); } Rprintf("\n");
		// for (j=0;j<dim_GARCH1;j++) { Rprintf("%lf\t",PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j]); } Rprintf("\n");

		logLikelihood(C,e,Nsize,np,orders,dims,
			CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
			Y_series,
			CURVALUE_m_series,CURVALUE_e_series,CURVALUE_h2_series,
			&CURVALUE_logLikelihood_s);

		inv(PRIOR_COV,&dim_GARCH_b1,PRIOR_inv_COV);
		for (k=0;k<dim_GARCH_b1;k++) {
			XY[k] = 0.0;
			for (l=0;l<dim_GARCH_b1;l++) {
				XY[k] += (PRIOR_inv_COV[k*dim_GARCH_b1+l]*PRIOR_AVE[l]);
				XX[k*dim_GARCH_b1+l] = PRIOR_inv_COV[k*dim_GARCH_b1+l];
				}
			}

		// for (k=0;k<dim_GARCH_b1;k++) { Rprintf("%lf\t",XY[k]); } Rprintf("\n");
		// for (k=0;k<dim_GARCH_b2;k++) { Rprintf("%lf\t",XX[k]); } Rprintf("\n");
			
		for (j=0;j<order_max;j++) {
			for (k=0;k<dim_GARCH_b1;k++) {
				Zi[j*dim_GARCH_b1+k] = 0.0;
				}
			}

		for (j=0;j<(*Nsize);j++) {
			e2[j+order_max] = CURVALUE_e_series[j+order_max] * CURVALUE_e_series[j+order_max];
			h4[j+order_max] = CURVALUE_h2_series[j+order_max] * CURVALUE_h2_series[j+order_max];
			index_c = C[j] - 1;
			for (k=0;k<dim_GARCH_b1;k++) {
				Zi[(j+order_max)*dim_GARCH_b1+k] = CURVALUE_h2_series[j+order_max-1-k];
				for (l=0;l<dim_GARCH_b1;l++) {
					Zi[(j+order_max)*dim_GARCH_b1+k] +=
						(CURVALUE_PARA_GARCHs[index_c*dim_GARCH1+dim_GARCH_a1+l]
						*Zi[(j+order_max-1-l)*dim_GARCH_b1+k]);
					}
				}
			if (index_c==i) {
				Yi[0] = e2[j+order_max] - CURVALUE_h2_series[j+order_max];
				for (k=0;k<dim_GARCH_b1;k++) {
					Yi[0] += (Zi[(j+order_max)*dim_GARCH_b1+k]
						*CURVALUE_PARA_GARCHs[index_c*dim_GARCH1+dim_GARCH_a1+k]);
					Xi[k] = Zi[(j+order_max)*dim_GARCH_b1+k];
					}
				for (k=0;k<dim_GARCH_b1;k++) {
					XY[k] += (Xi[k]*Yi[0]*0.5/h4[j+order_max]);
					for (l=0;l<dim_GARCH_b1;l++) {
						XX[k*dim_GARCH_b1+l] += (Xi[k]*Xi[l]*0.5/h4[j+order_max]);
						}
					}
				}
			}

		// for (k=0;k<dim_GARCH_b1;k++) { Rprintf("%lf\t",CURVALUE_PARA_GARCHs[i*dim_GARCH1+dim_GARCH_a1+k]); } Rprintf("\n");
			
		inv(XX,&dim_GARCH_b1,inv_XX);
		for (k=0;k<dim_GARCH_b1;k++) {
			Local_AVE[k] = 0.0;
			for (l=0;l<dim_GARCH_b1;l++) {
				Local_AVE[k] += (inv_XX[k*dim_GARCH_b1+l]*XY[l]);
				Local_COV[k*dim_GARCH_b1+l] = 10.0 * inv_XX[k*dim_GARCH_b1+l];
				// if (k==l) Local_COV[k*dim_GARCH_a1+l] = 0.1; else Local_COV[k*dim_GARCH_a1+l] = 0.0;
				}
			}
/*
		for (k=0;k<dim_GARCH_b1;k++) {
			Rprintf("%lf\t",XY[k]);
			Rprintf("| ");
			for (l=0;l<dim_GARCH_b1;l++) {
				Rprintf("%lf\t",XX[k*dim_GARCH_b1+l]);
				}
			Rprintf("| ");
			for (l=0;l<dim_GARCH_b1;l++) {
				Rprintf("%lf\t",inv_XX[k*dim_GARCH_b1+l]);
				}
			Rprintf("| ");
			Rprintf("%lf\t",Local_AVE[k]);
			Rprintf("| ");
			for (l=0;l<dim_GARCH_b1;l++) {
				Rprintf("%lf\t",Local_COV[k*dim_GARCH_b1+l]);
				}
			Rprintf("\n");
			}
*/
		ln_d_mv_norm(CURVALUE_PARA_GARCHs+i*dim_GARCH1+dim_GARCH_a1,Local_AVE,Local_COV,&dim_GARCH_b1,&CURVALUE_logProposal_s);
		ln_d_mv_norm(CURVALUE_PARA_GARCHs+i*dim_GARCH1+dim_GARCH_a1,PRIOR_AVE,PRIOR_COV,&dim_GARCH_b1,&CURVALUE_logPrior_s);
		kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s - CURVALUE_logProposal_s;

		kappa_MAXIMUM = 0.0;
		first_cal = 1;
		index = 0;
		positive = 0;
		count = 0;

		do {

			do {
				r_mv_norm(Local_AVE,Local_COV,&dim_GARCH_b1,PROPOSAL_PARA_GARCHs+i*dim_GARCH1+dim_GARCH_a1);
				positive = 1;
				for (j=0;j<dim_GARCH_b1;j++) {
					if (PROPOSAL_PARA_GARCHs[i*dim_GARCH1+dim_GARCH_a1+j]<0.0) { positive = 0; }
					}
				} while(positive==0);

			logLikelihood(C,e,Nsize,np,orders,dims,
				CURVALUE_PARA_ARs,PROPOSAL_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);
			ln_d_mv_norm(PROPOSAL_PARA_GARCHs+i*dim_GARCH1+dim_GARCH_a1,Local_AVE,Local_COV,&dim_GARCH_b1,&PROPOSAL_logProposal_s);
			ln_d_mv_norm(PROPOSAL_PARA_GARCHs+i*dim_GARCH1+dim_GARCH_a1,PRIOR_AVE,PRIOR_COV,&dim_GARCH_b1,&PROPOSAL_logPrior_s);
			kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s - PROPOSAL_logProposal_s;

/*
			for (k=0;k<dim_GARCH_b1;k++) { Rprintf("%lf\t",PROPOSAL_PARA_GARCHs[i*dim_GARCH1+dim_GARCH_a1+k]); }
			for (k=0;k<dim_GARCH_b1;k++) { Rprintf("%lf\t",CURVALUE_PARA_GARCHs[i*dim_GARCH1+dim_GARCH_a1+k]); }

			Rprintf("%lf\t",kappa_PROPOSAL);
			Rprintf("%lf\t",kappa_CURVALUE);
			Rprintf("%lf\t",kappa_PROPOSAL-kappa_CURVALUE);
			Rprintf("\n");
			Rprintf("\n");
*/

			if (first_cal==1) {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else {
					acc_prob = exp(kappa_PROPOSAL-kappa_CURVALUE);
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					first_cal = 0;
					}
				} else {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else if (kappa_PROPOSAL<kappa_MAXIMUM) {
					acc_prob = 0.0;
					index = 0;
					} else {
					acc_prob = 1.0
						- (1.0-exp(kappa_PROPOSAL-kappa_CURVALUE))
						/ (1.0-exp(kappa_MAXIMUM-kappa_CURVALUE));
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					}
				}
			
			if (index==1) {
				for (j=0;j<dim_GARCH1;j++) {
					CURVALUE_PARA_GARCHs[i*dim_GARCH1+j] =
					PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j];
					}
				} else {
				for (j=0;j<dim_GARCH1;j++) {
					PROPOSAL_PARA_GARCHs[i*dim_GARCH1+j] =
					CURVALUE_PARA_GARCHs[i*dim_GARCH1+j];
					}
				}

			count++; if (count>*max_R) { index = 1; }
			} while(index==0);

		}

	free(e2_tilde);
	free(i2_tilde);
	free(e2);
	free(h4);
	free(XX);
	free(XY);
	free(Yi);
	free(Xi);
	free(Zi);
	free(inv_XX);
	free(Local_AVE);
	free(Local_COV);
	free(PRIOR_AVE);
	free(PRIOR_COV);
	free(PRIOR_inv_COV);
	}

// garch parameters sampling step -------------- End

// ar parameters sampling step ----------------- Begin

void sample_AR_AMDR_MCMC(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *CURVALUE_PARA_ARs,		// [7]
	double *PROPOSAL_PARA_ARs,		// [8]
	double *CURVALUE_PARA_GARCHs,	// [9]
	double *PROPOSAL_PARA_GARCHs,	// [10]
	double *PRIOR_ARs,				// [11]
	double *PRIOR_GARCHs,			// [12]
	double *Y_series,				// [13]
	double *CURVALUE_m_series,		// [14]
	double *PROPOSAL_m_series,		// [15]
	double *CURVALUE_e_series,		// [16]
	double *PROPOSAL_e_series,		// [17]
	double *CURVALUE_h2_series,		// [18]
	double *PROPOSAL_h2_series,		// [19]
	int *Option,					// [20]
	double *AMDR_AVE_ARs,			// [21]
	double *AMDR_COV_ARs,			// [22]
	double *AMDR_AVE_GARCHs,		// [23]
	double *AMDR_COV_GARCHs,		// [24]
	int *Iter,						// [25]
	int *max_R						// [26]
	) {
	int i,j,k,l;
	double CURVALUE_logLikelihood_s,CURVALUE_logPrior_s,CURVALUE_logProposal_s;
	double PROPOSAL_logLikelihood_s,PROPOSAL_logPrior_s,PROPOSAL_logProposal_s;
	double kappa_CURVALUE,kappa_PROPOSAL,kappa_MAXIMUM;
	double acc_prob;
	double r_unif_s;
	int dim_AR1 = dims[0];
	int dim_AR2 = dims[1];
	int first_cal,index,count;
	
	double *Local_AVE; Local_AVE = (double *) calloc(dim_AR1,sizeof(double));
	double *Local_COV; Local_COV = (double *) calloc(dim_AR2,sizeof(double));
	double *PRIOR_AVE; PRIOR_AVE = (double *) calloc(dim_AR1,sizeof(double));
	double *PRIOR_COV; PRIOR_COV = (double *) calloc(dim_AR2,sizeof(double));
	
	for (i=0;i<(*np);i++) {
		for (j=0;j<dim_AR1;j++) {
			PROPOSAL_PARA_ARs[i*dim_AR1+j] =
			CURVALUE_PARA_ARs[i*dim_AR1+j];
			}
		}

	for (k=0;k<dim_AR1;k++) { PRIOR_AVE[k] = PRIOR_ARs[2*k]; }
	for (k=0;k<dim_AR2;k++) { PRIOR_COV[k] = 0.0; }
	for (k=0;k<dim_AR1;k++) { PRIOR_COV[k*dim_AR1+k] = PRIOR_ARs[2*k+1]; }

	for (i=0;i<(*np);i++) {

		for (j=0;j<dim_AR1;j++) {
			PROPOSAL_PARA_ARs[i*dim_AR1+j] =
			CURVALUE_PARA_ARs[i*dim_AR1+j];
			}

		// for (k=0;k<dim_AR1;k++) { Local_AVE[k] = CURVALUE_PARA_ARs[i*dim_AR1+k]; }
		for (k=0;k<dim_AR1;k++) { Local_AVE[k] = AMDR_AVE_ARs[i*dim_AR1+k]; }
		for (k=0;k<dim_AR2;k++) { Local_COV[k] = 2.34 * AMDR_COV_ARs[i*dim_AR2+k]; }

		logLikelihood(C,e,Nsize,np,orders,dims,
			CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
			Y_series,
			CURVALUE_m_series,CURVALUE_e_series,CURVALUE_h2_series,
			&CURVALUE_logLikelihood_s);
		ln_d_mv_norm(CURVALUE_PARA_ARs+i*dim_AR1,Local_AVE,Local_COV,&dim_AR1,&CURVALUE_logProposal_s);
		ln_d_mv_norm(CURVALUE_PARA_ARs+i*dim_AR1,PRIOR_AVE,PRIOR_COV,&dim_AR1,&CURVALUE_logPrior_s);
		kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s - CURVALUE_logProposal_s;
		// kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s;

		kappa_MAXIMUM = 0.0;
		first_cal = 1;
		index = 0;
		count = 0;

		do {

			r_mv_norm(Local_AVE,Local_COV,&dim_AR1,PROPOSAL_PARA_ARs+i*dim_AR1);

			logLikelihood(C,e,Nsize,np,orders,dims,
				PROPOSAL_PARA_ARs,CURVALUE_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);
			ln_d_mv_norm(PROPOSAL_PARA_ARs+i*dim_AR1,Local_AVE,Local_COV,&dim_AR1,&PROPOSAL_logProposal_s);
			ln_d_mv_norm(PROPOSAL_PARA_ARs+i*dim_AR1,PRIOR_AVE,PRIOR_COV,&dim_AR1,&PROPOSAL_logPrior_s);
			kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s - PROPOSAL_logProposal_s;
			// kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s;

			if (first_cal==1) {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else {
					acc_prob = exp(kappa_PROPOSAL-kappa_CURVALUE);
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					first_cal = 0;
					}
				} else {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else if (kappa_PROPOSAL<kappa_MAXIMUM) {
					acc_prob = 0.0;
					index = 0;
					} else {
					acc_prob = 1.0
						- (1.0-exp(kappa_PROPOSAL-kappa_CURVALUE))
						/ (1.0-exp(kappa_MAXIMUM-kappa_CURVALUE));
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					}
				}
			
			if (index==1) {
				for (j=0;j<dim_AR1;j++) {
					CURVALUE_PARA_ARs[i*dim_AR1+j] =
					PROPOSAL_PARA_ARs[i*dim_AR1+j];
					}
				} else {
				for (j=0;j<dim_AR1;j++) {
					PROPOSAL_PARA_ARs[i*dim_AR1+j] =
					CURVALUE_PARA_ARs[i*dim_AR1+j];
					}
				}

			count++; if (count>*max_R) { index = 1; }
			} while(index==0);

		for (k=0;k<dim_AR1;k++) {
			AMDR_AVE_ARs[i*dim_AR1+k]
				= AMDR_AVE_ARs[i*dim_AR1+k]
				+ (CURVALUE_PARA_ARs[i*dim_AR1+k]-AMDR_AVE_ARs[i*dim_AR1+k])
				/ (*Iter+1);
			}

		for (k=0;k<dim_AR1;k++) {
			for (l=0;l<dim_AR1;l++) {
				AMDR_COV_ARs[i*dim_AR2+k*dim_AR1+l]
					= AMDR_COV_ARs[i*dim_AR2+k*dim_AR1+l] +
					((CURVALUE_PARA_ARs[i*dim_AR1+k]-AMDR_AVE_ARs[i*dim_AR1+k])
					*(CURVALUE_PARA_ARs[i*dim_AR1+l]-AMDR_AVE_ARs[i*dim_AR1+l])
					-AMDR_COV_ARs[i*dim_AR2+k*dim_AR1+l]) / (*Iter+1);
				}
			}
		}

	free(Local_AVE);
	free(Local_COV);
	free(PRIOR_AVE);
	free(PRIOR_COV);
	}

void sample_AR_p_MCMC(
	int *C, 						// [1]
	int *e, 						// [2]
	int *Nsize, 					// [3]
	int *np, 						// [4]
	int *orders,					// [5]
	int *dims,						// [6]  // (dim.ar1,dim.ar2,dim.garch1,dim.garch2,dim.garch.a1,dim.garch.a2,dim.garch.b1,dim.garch.b2)
	double *CURVALUE_PARA_ARs,		// [7]
	double *PROPOSAL_PARA_ARs,		// [8]
	double *CURVALUE_PARA_GARCHs,	// [9]
	double *PROPOSAL_PARA_GARCHs,	// [10]
	double *PRIOR_ARs,				// [11]
	double *PRIOR_GARCHs,			// [12]
	double *Y_series,				// [13]
	double *CURVALUE_m_series,		// [14]
	double *PROPOSAL_m_series,		// [15]
	double *CURVALUE_e_series,		// [16]
	double *PROPOSAL_e_series,		// [17]
	double *CURVALUE_h2_series,		// [18]
	double *PROPOSAL_h2_series,		// [19]
	int *Option,					// [20]
	double *AMDR_AVE_ARs,			// [21]
	double *AMDR_COV_ARs,			// [22]
	double *AMDR_AVE_GARCHs,		// [23]
	double *AMDR_COV_GARCHs,		// [24]
	int *Iter,						// [25]
	int *max_R						// [26]
	) {
	int i,j,k,l;
	int index_c;
	double CURVALUE_logLikelihood_s,CURVALUE_logPrior_s,CURVALUE_logProposal_s;
	double PROPOSAL_logLikelihood_s,PROPOSAL_logPrior_s,PROPOSAL_logProposal_s;
	double kappa_CURVALUE,kappa_PROPOSAL,kappa_MAXIMUM;
	double acc_prob;
	double r_unif_s;
	int order_max = orders[0];
	int dim_AR1 = dims[0];
	int dim_AR2 = dims[1];
	int first_cal,index,count;

	double *h2; h2 = calloc((*Nsize)+order_max,sizeof(double));
	double *XY; XY = calloc(dim_AR1,sizeof(double));
	double *XX; XX = calloc(dim_AR2,sizeof(double));
	double *Yi; Yi = calloc(1,sizeof(double));
	double *Xi; Xi = calloc(dim_AR1,sizeof(double));
	double *inv_XX; inv_XX = calloc(dim_AR2,sizeof(double));
	double *Local_AVE; Local_AVE = (double *) calloc(dim_AR1,sizeof(double));
	double *Local_COV; Local_COV = (double *) calloc(dim_AR2,sizeof(double));
	double *PRIOR_AVE; PRIOR_AVE = (double *) calloc(dim_AR1,sizeof(double));
	double *PRIOR_COV; PRIOR_COV = (double *) calloc(dim_AR2,sizeof(double));
	double *PRIOR_inv_COV; PRIOR_inv_COV = (double *) calloc(dim_AR2,sizeof(double));

	for (i=0;i<(*np);i++) {
		for (j=0;j<dim_AR1;j++) {
			PROPOSAL_PARA_ARs[i*dim_AR1+j] =
			CURVALUE_PARA_ARs[i*dim_AR1+j];
			}
		}

	for (k=0;k<dim_AR1;k++) { PRIOR_AVE[k] = PRIOR_ARs[2*k]; }
	for (k=0;k<dim_AR2;k++) { PRIOR_COV[k] = 0.0; }
	for (k=0;k<dim_AR1;k++) { PRIOR_COV[k*dim_AR1+k] = PRIOR_ARs[2*k+1]; }

	for (i=0;i<(*np);i++) {

		for (j=0;j<dim_AR1;j++) {
			PROPOSAL_PARA_ARs[i*dim_AR1+j] =
			CURVALUE_PARA_ARs[i*dim_AR1+j];
			}

		logLikelihood(C,e,Nsize,np,orders,dims,
			CURVALUE_PARA_ARs,CURVALUE_PARA_GARCHs,
			Y_series,
			CURVALUE_m_series,CURVALUE_e_series,CURVALUE_h2_series,
			&CURVALUE_logLikelihood_s);

		inv(PRIOR_COV,&dim_AR1,PRIOR_inv_COV);
		for (k=0;k<dim_AR1;k++) {
			XY[k] = 0.0;
			for (l=0;l<dim_AR1;l++) {
				XY[k] += (PRIOR_inv_COV[k*dim_AR1+l]*PRIOR_AVE[l]);
				XX[k*dim_AR1+l] = PRIOR_inv_COV[k*dim_AR1+l];
				}
			}

		for (j=0;j<(*Nsize);j++) {
			index_c = C[j] - 1;
			if (index_c==i) {
				h2[j+order_max] = CURVALUE_h2_series[j+order_max];
				Yi[0] = Y_series[j+order_max];
				Xi[0] = 1.0;
				for (k=1;k<dim_AR1;k++) { Xi[k] = Y_series[j+order_max-k]; }
				for (k=0;k<dim_AR1;k++) {
					XY[k] += (Xi[k]*Yi[0]/h2[j+order_max]);
					for (l=0;l<dim_AR1;l++) {
						XX[k*dim_AR1+l] += (Xi[k]*Xi[l]/h2[j+order_max]);
						}
					}
				}
			}

		inv(XX,&dim_AR1,inv_XX);
		for (k=0;k<dim_AR1;k++) {
			Local_AVE[k] = 0.0;
			for (l=0;l<dim_AR1;l++) {
				Local_AVE[k] += (inv_XX[k*dim_AR1+l]*XY[l]);
				Local_COV[k*dim_AR1+l] = 10.0 * inv_XX[k*dim_AR1+l];
				}
			}

		ln_d_mv_norm(CURVALUE_PARA_ARs+i*dim_AR1,Local_AVE,Local_COV,&dim_AR1,&CURVALUE_logProposal_s);
		ln_d_mv_norm(CURVALUE_PARA_ARs+i*dim_AR1,PRIOR_AVE,PRIOR_COV,&dim_AR1,&CURVALUE_logPrior_s);
		kappa_CURVALUE = CURVALUE_logLikelihood_s + CURVALUE_logPrior_s - CURVALUE_logProposal_s;

		kappa_MAXIMUM = 0.0;
		first_cal = 1;
		index = 0;
		count = 0;

		do {

			r_mv_norm(Local_AVE,Local_COV,&dim_AR1,PROPOSAL_PARA_ARs+i*dim_AR1);

			logLikelihood(C,e,Nsize,np,orders,dims,
				PROPOSAL_PARA_ARs,CURVALUE_PARA_GARCHs,
				Y_series,
				PROPOSAL_m_series,PROPOSAL_e_series,PROPOSAL_h2_series,
				&PROPOSAL_logLikelihood_s);
			ln_d_mv_norm(PROPOSAL_PARA_ARs+i*dim_AR1,Local_AVE,Local_COV,&dim_AR1,&PROPOSAL_logProposal_s);
			ln_d_mv_norm(PROPOSAL_PARA_ARs+i*dim_AR1,PRIOR_AVE,PRIOR_COV,&dim_AR1,&PROPOSAL_logPrior_s);
			kappa_PROPOSAL = PROPOSAL_logLikelihood_s + PROPOSAL_logPrior_s - PROPOSAL_logProposal_s;

			if (first_cal==1) {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else {
					acc_prob = exp(kappa_PROPOSAL-kappa_CURVALUE);
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					first_cal = 0;
					}
				} else {
				if (kappa_CURVALUE<kappa_PROPOSAL) {
					acc_prob = 1.0;
					index = 1;
					} else if (kappa_PROPOSAL<kappa_MAXIMUM) {
					acc_prob = 0.0;
					index = 0;
					} else {
					acc_prob = 1.0
						- (1.0-exp(kappa_PROPOSAL-kappa_CURVALUE))
						/ (1.0-exp(kappa_MAXIMUM-kappa_CURVALUE));
					r_unif(&r_unif_s);
					if (r_unif_s<acc_prob) { index = 1; } else { index = 0; }
					kappa_MAXIMUM = kappa_PROPOSAL;
					}
				}
			
			if (index==1) {
				for (j=0;j<dim_AR1;j++) {
					CURVALUE_PARA_ARs[i*dim_AR1+j] =
					PROPOSAL_PARA_ARs[i*dim_AR1+j];
					}
				} else {
				for (j=0;j<dim_AR1;j++) {
					PROPOSAL_PARA_ARs[i*dim_AR1+j] =
					CURVALUE_PARA_ARs[i*dim_AR1+j];
					}
				}

			count++; if (count>*max_R) { index = 1; }
			} while(index==0);

		}

	free(h2);
	free(XY);
	free(XX);
	free(Yi);
	free(Xi);
	free(inv_XX);
	free(Local_AVE);
	free(Local_COV);
	free(PRIOR_AVE);
	free(PRIOR_COV);
	free(PRIOR_inv_COV);
	}

// ar parameters sampling step ----------------- End

