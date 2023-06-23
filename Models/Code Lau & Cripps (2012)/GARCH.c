#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "GARCH.h"

void sim_main(
	int *n, 			// [1]
	int *options, 		// [2]		
	int *nruns,			// [3]
	int *all_nps,		// [4]
	double *PARAs		// [5]
	) {
/*
	int i,j,k;

	int *Nsize; Nsize = (int *) calloc(1,sizeof(int)); *Nsize = *n;
	int process = *options; // 1 = PD; 2 = NGG;
	int nwarmup = *(nruns);
	int niter = *(nruns+1);
	int ntotalruns = nwarmup + niter;
	int *C; C = (int *) calloc(*Nsize,sizeof(int));
	int *e; e = (int *) calloc(*Nsize,sizeof(int));
	int *np; np = (int *) calloc(1,sizeof(int));
	
	// Rprintf("%d\t%d\t%d\n",process,sampling,steps);
	// Rprintf("%d\t%d\t%d\n",nreplication,nwarmup,niter);
	
	void (*sample_partition)(int *,int *,int *,int *,double *,double *,double *,double *,double *);

	if (process==1) {
		sample_partition = sample_PD_GIB_1STEP_MCMC;
		} else {
		sample_partition = sample_NGG_GIB_2STEP_MCMC;
		}
*/
/*
	for (j=0;j<nreplication;j++) {
		initial_partition(C,e,Nsize,np);
		for (i=0;i<ntotalruns;i++) {
			sample_partition(C,e,Nsize,np,PARAs);
			k = i - nwarmup;
			if (k>=0) {
				all_nps[j*(niter)+k] = *np;
				}
			}
		}
	free(C);
	free(e);
	free(np);
	free(Nsize);
*/
	}






