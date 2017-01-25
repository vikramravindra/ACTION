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

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iomanip>
using namespace std;

#include "CoinedPValues.h"

int main(int argc,char ** argv)
{
  int i,n,flag,flag_w;
  int N=1000;

  //random number generador seed
  int seed;

  //p-value cutoff
  double tau;
  double pvalue;
  
  //Vector used to store P-values and its weights
  struct COINED_PVALUE temp;
  vector<COINED_PVALUE> pvalue1;
  

  printf("\nCode description:\nThe current code computes the combine weighted P-value using (Bhoj,Fisher,Good, Mathai, Alves_Yu (cluster), Alves_Yu, RTGD, Gamma distribution) and compare it against the combine weighted P-value  obtained using the controlled approximation provided by Alves and Yi-Kuo (ArXiv e-print 1011.6627, 2010.<http://http://arxiv.org/abs/1011.6627>\n\n");   

   if(argc != 6)
      {

        printf("Combine_pvalue: <Number_of_Pvalues_to_combine> <Pvalue_weights:-1=Equal P-values and weights=1.0,0=Equal P-values differen weights,1==Equal weights,2==Unequal weights,3==Equal and unequal weights, 4==Equal weights>> <(0==Bhoj),(1==Fisher),(2==Good),(3==Mathai),(4==Alves & Yu (cluster)),(5==Alves & Yu),(6==Ratio of two Gamma Distribution (RTGD)), (7==Gamma Distribution), (8==TPM_Fisher)> <tau <=1> <seed>\n");
         exit(1);
      }

 n=atoi(argv[1]);
 flag_w=atoi(argv[2]);
 flag=atoi(argv[3]);
 tau=atof(argv[4]); //p-value cufoff value (truncate p-value cutoff)
 seed=atoi(argv[5]);

 srand48(seed);

 if(n <= 4) { n=4;}

 
  if(flag_w==-1) //Equal P-values equal weights
   {
      temp.pvalue=1.0/N;
      temp.pvalue_weight=1.0;
     for(i=0;i<n;i++)
       {
	 pvalue1.push_back(temp);
	  cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
       }
   }
 else if(flag_w==0) //Equal P-values different weights
   {
     temp.pvalue=1.0/N;
     for(i=0;i<n;i++)
       {
	 temp.pvalue_weight=drand48();
	 pvalue1.push_back(temp);
	  cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
       }
   }
 else if(flag_w==1) //Weights==1
   {
     for(i=0;i<n;i++)
       {
	 temp.pvalue=drand48()/N;
	 temp.pvalue_weight=1.0;
	 pvalue1.push_back(temp);
	  cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
       }
   }
 else if(flag_w==2 )  //Unequal weights
   {
     for(i=0;i<n;i++)
       {   
	 temp.pvalue=drand48()/N; 
	 temp.pvalue_weight=drand48();
	 pvalue1.push_back(temp);
	 cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
       }
   }
 else if(flag_w==3) //Equal and unequal weights
   {
     int k=0;
     k=int(0.40*n);  
     for(i=0;i<n-k;i++)
       {
	 temp.pvalue=drand48()/N;
	 temp.pvalue_weight=drand48();
	 pvalue1.push_back(temp);
	 cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
       }

      for(i=0;i<k;i=i+2)
       {
	 temp.pvalue  = drand48()/N;
	 pvalue1.push_back(temp);
	 cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
	 temp.pvalue  = drand48()/N;
	 pvalue1.push_back(temp);
	 cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
	 temp.pvalue_weight=drand48();
       }
   }
  else if(flag_w==4) //Equal 
   {
      temp.pvalue_weight=drand48();
     for(i=0;i<n;i++)
       {
	 temp.pvalue=drand48()/N;
	 pvalue1.push_back(temp);
	  cout<<"P-value = "<< temp.pvalue<<"\t"<<"Weight = "<<temp.pvalue_weight<<endl;
       }
   }

 CoinedPValues TEST(tau);   //creating an object CoinedPvalues
  
    if(flag==0)
      {    
        pvalue = TEST.CoinedPvalues_Log_Bhoj(pvalue1);      
    	cout<<"Bhoj's log(P-value) = "<<pvalue<<"\t P-value = "<<exp(pvalue)<<endl;
      }
    else if(flag==1)
      {
	pvalue = TEST.CoinedPvalues_Log_Fisher(pvalue1);
        cout<<"Fisher's log(P-value) = "<<pvalue<<"\t P-value = "<<exp(pvalue)<<endl;
      }
    else if(flag==2)
      {
        pvalue = TEST.CoinedPvalues_Log_Good(pvalue1);
        cout<<"Good's log(P-value) = "<<pvalue<<"\t P-value = "<<exp(pvalue)<<endl;
      } 
    else if(flag==3)
      {
        pvalue =  TEST.CoinedPvalues_Log_Mathai(pvalue1);
       cout<<"Mathai's log(P-value) = "<<pvalue<<"\t P-value="<<exp(pvalue)<<endl;
      } 
    else if(flag==4)
      {     
	pvalue = TEST.CoinedPvalues_Alves_Yu_Cluster(pvalue1); //Alves && Yu handles equal and unequal weights
	cout<<"Alves & Yu (cluster version) implemented P-value = "<<pvalue<<endl;
      }
    else if(flag==5)
      {     
	pvalue = TEST.CoinedPvalues_Alves_Yu(pvalue1); //Alves && Yu handles equal and unequal weights
	cout<<"Alves & Yu implemented P-value = "<<pvalue<<endl;
      }
    else if(flag==6)
      {  
	pvalue = TEST.CoinedPvalues_Log_WP_RTGD(pvalue1);
	cout<<"Ration of two Gamma Distribution P-value = "<<pvalue<<"\t P-value="<<exp(pvalue)<<endl; 
      }   
    else if(flag==7)
      {  
	pvalue = TEST.CoinedPvalues_Log_WP_GD(pvalue1);
	cout<<"Gamma Distribution log(P-value) = "<<pvalue<<"\t P-value="<<exp(pvalue)<<endl;
      }   
     else if(flag==8)
      {  
	pvalue = TEST.CoinedPvalues_Log_TPM_EFFD_Fisher(pvalue1);
	cout<<"TPM_Fisher Distribution log(P-value) = "<<pvalue<<"\t P-value="<<exp(pvalue)<<endl;
      }   
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
