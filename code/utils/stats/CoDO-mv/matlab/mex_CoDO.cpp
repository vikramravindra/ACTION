#include <math.h>
#include "mex.h"

void CoDO(int x, int nL, int *L, int n_union, int n, int m_overlap, int m_union, double *p);

int min(int x, int y) {
	return (x < y? x:y);
}

int max(int x, int y) {
	return (x > y? x:y);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	register int i, j;
    int row, col;
	
	int NumThreads = 1;
	if(nrhs < 6) {
		mexErrMsgTxt("Syntax: mex_CoDO(set_sizes, overlap_size, union_size, population_size, overlap_edges, union_edges)\n");					
		
	}
	
	int verbose = 0;
	if(nrhs == 7)
		verbose = (int ) mxGetScalar(prhs[6]);
		

	if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgTxt("List of set sizes should be double\n");
	}

	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);	
	
	int min_sizes = min(m, n);
	int max_sizes = max(m, n);
	
	if(min_sizes != 1 || max_sizes < 2) {
		mexErrMsgTxt("List of set sizes should be a vector with size at least 2\n");		
	}	

	if(verbose)
		mexPrintf("CoDO with %d sets\n", max_sizes);

	double *L_dbl = (double *) mxGetPr(prhs[0]);
	int *L = new int[max_sizes];
	for (i = 0; i < max_sizes; i++) {
		L[i] = (int)roundl(L_dbl[i]);
	}
	
	for (i = 0; i < max_sizes; i++) {
		if(verbose)		
			mexPrintf("\tSet %d size = %d\n", i, L[i]);
	}

	
	if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
		mexErrMsgTxt("overlap_size should be an integer\n");					
	}
	int overlap_size = (int ) mxGetScalar(prhs[1]);
	if(verbose)
		mexPrintf("Overlap size = %d\n", overlap_size);
		
	
	if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
		mexErrMsgTxt("union size should be an integer\n");					
	}
	int union_size = (int ) mxGetScalar(prhs[2]);
	if(verbose)
		mexPrintf("union size = %d\n", union_size);



	if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ) {
		mexErrMsgTxt("population size should be an integer\n");					
	}
	int pop_size = (int ) mxGetScalar(prhs[3]);
	if(verbose)
		mexPrintf("population size = %d\n", pop_size);
		

	if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ) {
		mexErrMsgTxt("overlap edges should be an integer\n");					
	}
	int overlap_edges = (int ) mxGetScalar(prhs[4]);
	if(verbose)
		mexPrintf("overlap edges = %d\n", overlap_edges);


	if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ) {
		mexErrMsgTxt("union edges should be an integer\n");					
	}
	int union_edges = (int ) mxGetScalar(prhs[5]);
	if(verbose)
		mexPrintf("union edges = %d\n", union_edges);
		
	mexEvalString("drawnow;");	
	double pval;
	if(overlap_edges < 1) {
		mexPrintf("empty edge set -> pval = 1\n");
		pval = 1;
	}
	else
		CoDO(overlap_size, max_sizes, L, union_size, pop_size, overlap_edges, union_edges, &pval);
			
	if(pval > 1) pval = 0; // Overflow!
	
	if(verbose)
		mexPrintf("final pval = %e\n", pval);
	
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(plhs[0]) = pval;
}
