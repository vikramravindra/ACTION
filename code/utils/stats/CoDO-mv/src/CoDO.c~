#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mvhyper.h"

void CoDO(int x, int nL, int *L, int n_union, int n, int m_overlap, int m_union, double *p){
/*
x:         number of elements overlap between all subsets
nL:        number of subsets
L:         subset sizes
n_union:   size of union set
n:         background size
m_overlap: Number of edges in intersection
m_union:   Number of edges in union
p:         output probability
*/
	const double tiny = 1.0E-320;
	register int i;
	
	int minL = min(L,nL);		


	int total_pairs = (int) ( n_union*(n_union - 1) / 2.0 ) ;
	
	double Xmean = 0.0 + n;
	for(i = 0; i < nL ; i++){
		Xmean = Xmean * L[i] / n;
	}
		
	
	int densityL[2], n_densityL = 2;;
	densityL[0] = m_union; // sample size in hygecdf
	double density_pval=1;
	int success_size = m_overlap - 1;
	

	// Precompute and store log values
	int max_size = (total_pairs < n? n:total_pairs);	
	double* logVal = (double *)calloc(max_size, sizeof(double));
	if(logVal == NULL) {
		fprintf(stderr, "CoDO: logVal memory allocation failed\n"); fflush(stderr);
	}	
	for(i=1; i<= max_size ; i++){
		logVal[i-1]=log((double)i);
	}


	int lower = 0, logp = 0;
	double pval = 0.0, size_prob=0.0;


	double last_p = 1.0;
	for(i = x; i <= minL; i++){			
		C_dmvhyper_logVal(&i, &nL, L, &n, &size_prob, &logp, logVal);					
		
		if(size_prob <= tiny) break;
		if( size_prob/last_p <= 1e-3) break;  //No improve in precision
		last_p = size_prob;
		
		if(i <= 1) 
			density_pval = 1;
		else {
			densityL[1] = i*(i-1) / 2; // pos size in hygecdf
			
			// sample size = L[0] = m_union, pos size = L[1] = C(i, 2), success_size = m_overlap - 1
			C_pmvhyper_logVal(&success_size, &n_densityL, densityL, &total_pairs, &density_pval, &lower, &logp, logVal);							
		}
		pval += size_prob*density_pval;	
		printf(" size prob = %e, density pval = %d, pval = %e\n", size_prob, density_pval, pval);	
	}
	
	if(1 < pval) {
		*p = 1.0;
	}
	else {	
		(*p) = pval;
	}

	printf("final pval = %e\n", *p);
	
	
	free(logVal);
	return;
}
