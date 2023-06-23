#ifndef CCRPM_H_
#define CCRPM_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

void r_gamma(double *,double *,double *);
void ln_r_unif(double *);
void r_unif(double *);
void chol(double *,int *,double *);
void chol2inv(double *,int *,double *);
void inv(double *,int *,double *);
void eigen(double *,int *,double *,double *);
void ln_d_t(double *,double *,double *,double *,double *);
void ln_d_Ltnorm(double *,double *,double *,double *);
void ln_d_norm(double *,double *,double *,double *);
void ln_d_uv_norm(double *,double *,double *,double *);
void ln_d_mv_norm(double *,double *,double *,int *,double *);
void r_norm(double *,double *,double *);
void r_mv_norm_chol(double *,double *,int *,double *);
void r_mv_norm(double *,double *,int *,double *);
void r_uv_norm(double *,double *,int *,double *);
void r_Ltnorm(double *,double *,double *);
double nrs_a_inf(double *);
double ers_a_inf(double *);

void initial_GIB_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *);
void Pred_PD_GIB_1STEP_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,int *);
void Pred_PD_GIB_2STEP_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,int *);
void Pred_NGG_GIB_2STEP_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,int *);
void sample_PD_GIB_1STEP_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *);
void sample_PD_GIB_2STEP_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *);
void sample_NGG_GIB_2STEP_MN(int *,int *,int *,int *,double *,double *,double *,double *,double *);

void logLikelihood(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *);
void logLikelihood_i_to_n(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,int *);
void logLikelihood_pred_PD(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *);
void logLikelihood_pred_NGG(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *);
void logLikelihood_pred_NON(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *);

void sample_PD_GIB_1STEP_MGARCH(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,int *,int *);
void sample_NON_GIB_1STEP_MGARCH(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,int *,int *);
void sample_NGG_GIB_2STEP_MGARCH(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,int *,int *);

void random_GARCH_pq(int*,int*,double *,double *);
void random_AR_p(int*,int*,double *,double *);

void sample_GARCH_pq_MCMC(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,double *,double *,double *,double *,int *,int *);
void sample_GARCH_AMDR_MCMC(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,double *,double *,double *,double *,int *,int *);
void sample_AR_p_MCMC(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,double *,double *,double *,double *,int *,int *);
void sample_AR_AMDR_MCMC(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,double *,double *,double *,double *,int *,int *);

void logcumsum(double *,double *,int *);


#endif

