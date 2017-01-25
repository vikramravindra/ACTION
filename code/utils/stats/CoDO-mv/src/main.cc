#include <stdio.h>

void CoDO(int x, int nL, int *L, int nU, int n, int m_overlap, int m_union, double *p);
/* -- Test 1 --
 * Output:
	Size Set 0 = 38
	Size Set 1 = 31
	Overlap size = 30, Union size = 39, n = 2246
	Overlap edges = 175, union edges = 294
	pval = 5.388694e-60
 */

/* -- Test 2 --
 * Output:
	Size Set 0 = 118
	Size Set 1 = 249
	Overlap size = 117, Union size = 250, n = 2639
	Overlap edges = 2128, union edges = 8159
	pval = 4.181543e-157
 */

 
/* -- Test 3 --
 * Output:
	Size Set 0 = 645
	Size Set 1 = 473
	Size Set 2 = 309
	Size Set 3 = 79
	Overlap size = 8, Union size = 909, n = 4645
	Overlap edges = 28, union edges = 231103
	pval = 5.311639e-19
*/

int main() {
//  Test 1

/*	
	const int size = 2;
	const int nZ = 30, n = 2246, nU = 39;
	const int mZ = 175, mU = 294;
	int L[size] = {38, 31};
*/	
	
//  Test 2
/*
	const int size = 2;
	const int nZ = 117, n = 2639, nU = 250;
	const int mZ = 2128, mU = 8159;
	int L[size] = {118, 249};
*/


// Test 3
/*
	const int size = 4;
	const int nZ = 8, n = 4645, nU = 909;
	const int mZ = 28, mU = 231103;
	int L[size] = {645, 473, 309, 79};
*/	
	
	
	const int size = 2;
	const int nZ = 32, n = 4645, nU = 1037;
	const int mZ = 468, mU = 215184;
	int L[size] = {817, 252};
	
		
	double pval;
	CoDO(nZ, size, L, nU, n, mZ, mU, &pval);
	
	for(int i = 0; i < size; i++) {
		printf("Size Set %d = %d\n", i, L[i]);
	}
	printf("Overlap size = %d, Union size = %d, n = %d\n", nZ, nU, n);
	printf("Overlap edges = %d, union edges = %d\n", mZ, mU);
	printf("pval = %e\n", pval);
	
	return 0;	
}
