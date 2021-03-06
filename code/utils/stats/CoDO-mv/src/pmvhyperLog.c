#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mvhyper.h"


void C_pmvhyper_logVal(int *x, int *nL, int *L, int *n, double *p, int *lower, int *logp, double* logVal){
/*
x:     number of elements overlap between all subsets
nL:    number of subsets
L:     subset sizes
n:     background size
p:     output probability
lower: 1, lower tail probability Pr(overlap <= x); 0, upper tail probability Pr(overlap > x)
logp:  return log probability
*/
	const double tiny = 1.0E-320;
	int i,j;
	int i0=0;
	double p0=0.0;


	double Xmean;
	int minL=min(L,*nL);
	double* pp = (double *)calloc(minL+1, sizeof(double));


	if(*x == 0){
		C_dmvhyper_logVal(x, nL, L, n, p, &i0, logVal);
		if(*lower == 0) *p = 1.0 - *p;
		if(*logp>0) *p=log(*p);
		return;
	}
	Xmean=0.0 + *n;
	for(i=0; i< *nL ; i++){
		Xmean = Xmean * L[i] / *n;
	}
	*p = 0.0;
	if((double) *x > Xmean){
		i = *x + 1;
		for(; i <= minL; i++){
			C_dmvhyper_logVal(&i, nL, L, n, &p0, &i0, logVal);
			pp[i]=p0;
			if(p0 <= tiny) break;
			if(i > (*x + 1) && (p0/pp[i-1]) <= 0.01) break;  //No improve in precision
		}
		for(j = i; j >= *x + 1; j--){ //iteration from smallest to largest; more accurate
			*p += pp[j];
		}
		if(*lower > 0) *p = 1.0 - *p;
	}else{
		i = *x;
		for(;i >= 0; i--){
			C_dmvhyper_logVal(&i, nL, L, n, &p0, &i0, logVal);
			pp[i]=p0;
			if(p0 <= tiny) break;
			if(i < *x && (p0/pp[i+1]) <= 0.01) break;  //No improve in precision
		}
		for(j = i; j <= *x; j++){ //iteration from smallest to largest; more accurate
			*p += pp[j];
		}
		if(*lower == 0) *p = 1.0 - *p;
	}
	
	if(*logp > 0) *p = log(*p);

	free(pp);
	return;
}
