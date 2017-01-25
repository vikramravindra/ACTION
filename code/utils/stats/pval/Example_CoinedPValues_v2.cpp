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
  int flag;
  FILE *input;
  //random number generador seed

  //p-value cutoff
  double tau;
  double pvalue;
  
  //Vector used to store P-values and its weights
  struct COINED_PVALUE temp;
  vector<COINED_PVALUE> pvalue1;
  

  printf("\nCode description:\nThe current code computes the combine weighted P-value using (Bhoj,Fisher,Good, Mathai, Alves_Yu (cluster), Alves_Yu, RTGD, Gamma distribution) and compare it against the combine weighted P-value  obtained using the controlled approximation provided by Alves and Yi-Kuo (ArXiv e-print 1011.6627, 2010.<http://http://arxiv.org/abs/1011.6627>\n\n");   

   if(argc != 4)
      {

        printf("Combine_pvalue: <P-value File> <(0==Bhoj),(1==Fisher),(2==Good),(3==Mathai),(4==Alves & Yu (cluster)),(5==Alves & Yu),(6==Ratio of two Gamma Distribution (RTGD)), (7==Gamma Distribution), (8==TPM_Fisher)> <tau <=1>\n");
         exit(1);
      }

  input=fopen(argv[1],"r");
  flag=atoi(argv[2]);
  tau=atof(argv[3]); //p-value cufoff value (truncate p-value cutoff)

 
  if(input==NULL) 
    {  
      cout<<"Can't open file = "<<argv[1]<<endl; exit(1);
    }
  else
    {
      while(fscanf(input,"%s  %lf  %lf\n",temp.item,&temp.pvalue,&temp.pvalue_weight) != EOF)
	{
           pvalue1.push_back(temp);
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
	cout<<"Alves & Yu implemented P-value = "<<pvalue<<"\t"<<log(pvalue)<<endl;
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
