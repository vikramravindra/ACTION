/*************************************************************************************************************************************************/
/*
Author: LCDR Gelio Alves, PhD.  Member: United States Public Health Services.
Group: QMBP/NCBI/NLM/NIH
Code start date: December, 2010.
Code finish date: January,2011
Code last updated: May 9,2011
Description: This code is used to combine weighted P-values using a controlable approximation.
References: ArXiv e-print 1011.6627, 2010. <http://http://arxiv.org/abs/1011.6627>
*/
/************************************************************************************************************************************************/
#include <iostream>
#include<iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <limits.h>
#include <float.h>
#include <complex>
#include <algorithm>
using namespace std;

//Store users input list of P-values
typedef struct COINED_PVALUE
{
  int flag;
  double pvalue;  //store pvalues
  double pvalue_weight;  //store p-value weight
  int n_pvalue; //store the number of p-values in the cluster
  char item[128]; //item key(name)
  inline bool operator < (const struct COINED_PVALUE &a) const {   return pvalue_weight < a.pvalue_weight; };
} COINED_PVALUE;

//Sorting using the pvalue as the variable 
struct sort_cpvalue
{
  //sorting increase order
  inline bool operator()(const struct  COINED_PVALUE &A,const struct  COINED_PVALUE &B) const { return A.pvalue < B.pvalue;} 
};

class CoinedPValues
{
 
 public:
 
  //Parameters:(Truncate product P-value cutoff)
  CoinedPValues(double tau=1.0);
  
  //Returns the combined  P-value for a list of P-values with weights using our formula (Gelio Alves and Yi-Kuo Yu) 
  //Reference:ArXiv e-print 1011.6627, 2010. <http://http://arxiv.org/abs/1011.6627>  
  //This function computes the combine P-value by first clustering P-values of similar weights.
  //This is done to speed up the computation of the unified P-value
  double CoinedPvalues_Alves_Yu_Cluster(vector<struct COINED_PVALUE>);  //Cluster identical weights before combine P-values 

  //Original function published at: ArXiv e-print 1011.6627, 2010. <http://http://arxiv.org/abs/1011.6627> 
  //Gelio alves and Yi-Kuo. Yu. Combining independent, weighted P-values: achieving computational stability by a systematic expansion with controllable accuracy. PLoS One. 2011;6(8):e22647. doi: 10.1371/journal.pone.0022647. Epub 2011 Aug 31.
  double CoinedPvalues_Alves_Yu(vector<struct COINED_PVALUE>);   //Do not cluster identical weights before combine P-values   

  //Returns the combined P-value for a list of P-values with their associated Weights. Using Mathai's formula.  
  //Reference:Mathai, A. (1983). On linear combinations of independent exponential variables. Communications in Statistics - Theory and Methods 12, 625-632.
  double CoinedPvalues_Mathai(vector<struct COINED_PVALUE>);
  double CoinedPvalues_Log_Mathai(vector<struct COINED_PVALUE>);

 //Returns the combined P-value for a list of P-values with their associated Weightes. Using Good's formula
 //Reference: Good, I. J. (1955). On the weighted combination of significance tests. Journal of the Royal  Statistical Society. Series B (Methodological) 17(2), 264265
  double CoinedPvalues_Good(vector<struct COINED_PVALUE>); 
  double CoinedPvalues_Log_Good(vector<struct COINED_PVALUE>);   

  //Returns the combined P-value for a list of P-values with their associated Weightes. Using Fisher's formula
  //Reference:Bailey, T. L. and M. Gribskov (1998). Combining evidence using p-values: application to sequence homology searches. Bioinformatics 14, 4854.
  double CoinedPvalues_Fisher(vector<struct COINED_PVALUE>); 
  double CoinedPvalues_Log_Fisher(vector<struct COINED_PVALUE>); 
  //Returns the P-value for the weighted combine truncated P-value using the effective dimension equal to the sum of weights
  double CoinedPvalues_Log_TPM_EFFD_Fisher(vector<struct COINED_PVALUE>);
  //Returns the P-value for the weighted combine truncated P-value using the effective dimension equal to the sum of weights
  double CoinedPvalues_Log_TPM_EFFD_Fisher_Tau1(vector<struct COINED_PVALUE>);
  //Combine fisher P-value with tau=1
  double CoinedPvalues_Fisher_Tau1(vector<struct COINED_PVALUE>);  
  double CoinedPvalues_Fisher_Log_Tau1(vector<struct COINED_PVALUE>);  

  //Return the combine P-value using Bhoj approximation
  //Reference:Bhoj, Dinesh. Statistics and Probability Letters 15 (1992) 37-40. 
  double CoinedPvalues_Bhoj(vector<struct COINED_PVALUE>);
  double CoinedPvalues_Log_Bhoj(vector<struct COINED_PVALUE>);
 
  //Combined Weighted P-values (WP) using ratio of two gamma distributions RTGD
  double CoinedPvalues_WP_RTGD(vector<struct COINED_PVALUE>);
  double CoinedPvalues_Log_WP_RTGD(vector<struct COINED_PVALUE>);

  //Combined Weighted P-values (WP) using gamma distributions GD
  double CoinedPvalues_WP_GD(vector<struct COINED_PVALUE>);
  double CoinedPvalues_Log_WP_GD(vector<struct COINED_PVALUE>);

  //Combining TPM Zaykin
  //Reference: Dmitri V Zaykin and etal. Truncated Method for Combining p-value. Gent Epidemiol. 2002 Feb;22(2):170-85.
  double CoinedPvalues_TPM_Zaykin(vector<struct COINED_PVALUE>);

  //Needs to implement Moschopoulos
  //Reference: P. G. Moschopoulos (1985). The Distribution of The Sum of Independent Gamma Random Variables. 37 (1985), Part A, 541-544.
  //double CoinedPvalues_Moschopoulos(vector<struct COINED_PVALUE>);
 
 private:


  //Removing degenerate weights from data by combining P-value together.
  void Remove_Degenerate_Weights(vector<struct COINED_PVALUE> &);

  //Ensure users P-values are between [0,1]
  void Check_Pvalue_Range(vector<struct COINED_PVALUE> &);

  //Computes the P-value for the weighted combined P-value
  double Combine_Weighted_Pvalue_CDF(vector<struct CCW_COINED_PVALUE>,double,double);

  //Computes the F-terms  terms for the CDF 
  double Combine_Weighted_Pvalue_Approx(vector<struct CCW_COINED_PVALUE>,double); 

  //Recursion function for combine weight p-value CDF
  double Recursion_Combine_Weighted_Pvalue_CDF(int,int,int,int,vector<int>,vector<struct CCW_COINED_PVALUE>);  
  
  //Computes the coefficients Y-terms for correction terms
  double inline Combine_Weight_Pvalue_Delta_Weight(int,int,vector<struct CCW_COINED_PVALUE>);

  //Recursion function for combine weight p-value CDF
  double Recursion_Combine_Weighted_Pvalue_CDF_Mathai(int,int,int,int,vector<double>,vector<int> ,vector<int> );

  //Return the denominator term of Goods Formula
  void Good_Denominator(int, vector<struct COINED_PVALUE> &, long double &,int &);

  //Return P-value of Beta Distribution
  double Pvalue_Beta_Pdf(double a, double b, double x);
  double betacf(double a, double b, double x);

  //Used to compute incomplete gamma function
  void ln_gcf(double *, double, double, double *);
  void ln_gser(double *, double, double, double *);
  double ln_gammq(double,double); //Q(a,x) = 1- P(a,x)

  //Computes integer factorial
  void Factorial();

  //returns ln(n!)
  double FACTR_LN(int);
  //returns ln[gamma(xx)]  
  double GAMMA_LN(double xx);
  //Store values of factorial function
  double factorial[5001];


  double TAU;   //Truncate product P-value cutoff
};

//Store information used to numerically compute the combine P-value
 struct CCW_COINED_PVALUE
  {
    double pvalue;
    double pvalue_weight;
    int    pvalue_degeneracy;
    vector<double> pvalue_weight_dist;
    inline bool operator < (const struct CCW_COINED_PVALUE &a) const {   return pvalue_weight < a.pvalue_weight; };
 };
