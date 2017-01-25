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
#include "CoinedPValues.h"
CoinedPValues::CoinedPValues(double tau)
{
  TAU = tau;
  Factorial(); 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Removing degenerate weights from data.*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoinedPValues::Remove_Degenerate_Weights(vector<struct COINED_PVALUE> &v_pvalue)
{
  int i,j;
  double epsilon = 1e-6;      //epsilon to cluster identical weights
  double s_weight;
  vector<struct COINED_PVALUE>  c_pvalue;

  for(i=0;i<int(v_pvalue.size());i++)  
    {
      c_pvalue.clear(); 
      s_weight = 0;
      s_weight = v_pvalue[i].pvalue_weight;
      for(j=i+1;j<int(v_pvalue.size());j++)  
	{  	 
	  if( fabs(v_pvalue[i].pvalue_weight - v_pvalue[j].pvalue_weight) <= epsilon)
	    {
	      // s_weight = s_weight + v_pvalue[j].pvalue_weight; 
	      //v_pvalue[j].pvalue_weight = v_pvalue[i].pvalue_weight;
	      c_pvalue.push_back(v_pvalue[j]);
	      v_pvalue.erase(v_pvalue.begin() + j);
	      j = j - 1;
	    }  
	}
      if( c_pvalue.size() > 0 )
	{
	  c_pvalue.push_back(v_pvalue[i]);
	  //v_pvalue[i].pvalue = CoinedPvalues_Fisher(c_pvalue);
	  v_pvalue[i].pvalue = CoinedPvalues_Fisher_Tau1(c_pvalue);
	  v_pvalue[i].n_pvalue = c_pvalue.size();
	  v_pvalue[i].pvalue_weight = s_weight; 
	  //v_pvalue[i].pvalue = CoinedPvalues_Alves_Yu(c_pvalue);
          //When all the weights are equal Fisher P-value should be equal to Alves_Yu P-value
  //cout<<"CoinePvalue::Group p-value using Alves Yu = "<<CoinedPvalues_Alves_Yu(c_pvalue)<<", weight = "<<v_pvalue[i].pvalue_weight<<", group size = "<<c_pvalue.size()<<endl;
	  // cout<<"CoinedPvalue::Grouping P-value using Fisher="<<v_pvalue[i].pvalue<<", weight="<<v_pvalue[i].pvalue_weight<<", group size="<<c_pvalue.size()<<endl;
	}
    }
  //cout<<"CoinedPvalue::Vector size = "<<v_pvalue.size()<<endl;
  //for(i=0;i<int(v_pvalue.size());i++) {cout<<v_pvalue[i].pvalue<<"\t"<<v_pvalue[i].pvalue_weight<<endl;}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Returns the combined  P-value for a list of P-values with weights using our formula (Gelio Alves and Yi-Kuo Yu). 
Reference:ArXiv e-print 1011.6627, 2010. <http://http://arxiv.org/abs/1011.6627>*/ 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Alves_Yu_Cluster(vector<struct COINED_PVALUE> v_pvalue)
{
  double p;
 
  //Combining P-value with degenrate weights
  Remove_Degenerate_Weights(v_pvalue);

  if(  v_pvalue.size() > 1 )
    {
      p = CoinedPvalues_Alves_Yu(v_pvalue); 
    }
  else if( v_pvalue.size() == 1)
    {
      p = v_pvalue[0].pvalue; 
    }
  else
    {
      p=1;
    }

  return p;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Returns the combined  P-value for a list of P-values with weights using our formula (Gelio Alves and Yi-Kuo Yu). 
Reference:ArXiv e-print 1011.6627, 2010. <http://http://arxiv.org/abs/1011.6627>*/ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Alves_Yu(vector<struct COINED_PVALUE> v_pvalue)
{
  int i,j,k,flag;
  double t_weight;
  double Epsilon_Cluster=1e-6;  //to cluster close weights
  double epsilon = 1e-12;       //to cluster identical weights
  double Log_Tau=0;
  double Log_Coeff=0;
  double normalize_weights; 
  double final_pvalue=1.0;
  struct CCW_COINED_PVALUE temp;
  vector<struct CCW_COINED_PVALUE> pvalue;

  Check_Pvalue_Range(v_pvalue); //Checking P-value between [0,1]

////////Step 0:Initialize P-value structure and computes normalization constant
  //cout<<"Number of p-values = "<<v_pvalue.size()<<endl; 
  normalize_weights=0;
  for(i=0;i<int(v_pvalue.size());i++)  
    {   
      temp.pvalue =  v_pvalue[i].pvalue;
      temp.pvalue_weight = v_pvalue[i].pvalue_weight;
      temp.pvalue_degeneracy = 1;
      pvalue.push_back(temp);
      normalize_weights =  normalize_weights + 1.0/v_pvalue[i].pvalue_weight;
    }
  
/////////Step 1:Computes normalized P-value weightes and computes tau
  //cout<<"Number of p-values = "<<pvalue.size()<<endl; 
  Log_Tau=0;
  for(i=0;i<int(pvalue.size());i++)  
    {   
      pvalue[i].pvalue_weight = (1.0/pvalue[i].pvalue_weight)*(int(pvalue.size())/normalize_weights);			 
      pvalue[i].pvalue_weight_dist.push_back(pvalue[i].pvalue_weight);
      Log_Tau = Log_Tau - (1.0/pvalue[i].pvalue_weight)*log(pvalue[i].pvalue);
      // cout<<i<<"\t"<<pvalue[i].pvalue<<"\t"<<normalize_weights<<"\t"<<1.0/pvalue[i].pvalue_weight<<"\t"<<Log_Tau<<"\t"<<Log_Tau<<"\t"<<pow(pvalue[i].pvalue,1.0/pvalue[i].pvalue_weight)<<endl;
    }
  //cout<<"Step 1:Debug: Weighted P-value product: Log_Tau = "<<Log_Tau<<endl<<endl;

////////Step 2:Sort P-value 1/weightes in increasing order of weight
  
  sort(pvalue.begin(),pvalue.end());
  //cout<<"Step 2:Debug: sorting weight in increasing order"<<endl<<endl;
  
/////////Step 3:Cluster P-value 1/weightes in increasing order of weight
 
  //cout<<"Step 3:Debug: clustering P-values"<<endl;
  
  //Clustering identical P-values weights
  //cout<<"CoinedPValues::CoinedPvalues_Alves_Yu:Clustering identical p-values:"<<pvalue.size()<<endl;
  for(i=0;i<int(pvalue.size());i++)  
    {
      t_weight = pvalue[i].pvalue_weight;
      for(j=i+1;j<int(pvalue.size());j++)  
	{
	  if( fabs(t_weight - pvalue[j].pvalue_weight) <= epsilon)
	    {	      
	      pvalue[i].pvalue_weight = pvalue[i].pvalue_weight +  pvalue[j].pvalue_weight;
	      pvalue[i].pvalue_degeneracy= pvalue[i].pvalue_degeneracy + pvalue[j].pvalue_degeneracy; 
	      pvalue[i].pvalue_weight_dist.push_back(pvalue[i].pvalue_weight); 
	      pvalue.erase(pvalue.begin()+j,pvalue.begin()+j+1);	      
	      j=j-1;
	    }
	}
      pvalue[i].pvalue_weight = pvalue[i].pvalue_weight/pvalue[i].pvalue_degeneracy; 
    }
  //cout<<"CoinedPValues::CoinedPvalues_Alves_Yu:Clustering identical p-values:"<<pvalue.size()<<"\t"<<pvalue[i].pvalue_weight<<endl;
  //cout<<"Step 3.1:Debug: clustering P-values"<<endl;

  //Clustering P-values weights that are within Epsilon
  flag=1;
  while(flag == 1)
    {
      flag=0; 
      for(i=0;i<int(pvalue.size()-1);i++) 
	{
	  j=i+1;
           
	  if(j < int(pvalue.size()))
	    {
	      if( fabs(pvalue[i].pvalue_weight - pvalue[j].pvalue_weight) <=  Epsilon_Cluster)
		{
		  flag=1;

                  pvalue[i].pvalue_weight = ( pvalue[i].pvalue_weight*pvalue[i].pvalue_degeneracy  +  pvalue[j].pvalue_weight*pvalue[j].pvalue_degeneracy)/(pvalue[i].pvalue_degeneracy + pvalue[j].pvalue_degeneracy);
                  pvalue[i].pvalue_degeneracy = pvalue[i].pvalue_degeneracy + pvalue[j].pvalue_degeneracy; 

                  for(k=0;k<int(pvalue[j].pvalue_weight_dist.size());k++)
		    {
		      pvalue[i].pvalue_weight_dist.push_back(pvalue[j].pvalue_weight_dist[k]); 
		    }
		  // cout<<"i ="<<i<<"  "<<pvalue[i].pvalue_weight<<"\t"<<pvalue[i+1].pvalue_weight<<"\t"<<pvalue.size()<<"\t"<<"size = "<<pvalue[i].pvalue_weight_dist.size()<<"\t"<<"eta = "<<pvalue[i].pvalue_weight - pvalue[j].pvalue_weight<<endl;
		  pvalue.erase(pvalue.begin()+j,pvalue.begin()+j+1);
		}
	    }
	}
    }

/////////////Step 4:Computes eta
    
  for(i=0;i<int(pvalue.size());i++)
    {
      for(j=0;j<int(pvalue[i].pvalue_weight_dist.size());j++)
	{
	  // cout<<"i = "<<i<<"\t"<<j<<"\t"<< pvalue[i].pvalue_weight<<"\t"<<pvalue[i].pvalue_weight_dist[j]<<"\t"<<pvalue[i].pvalue_weight - pvalue[i].pvalue_weight_dist[j]<<endl;
	   pvalue[i].pvalue_weight_dist[j] = pvalue[i].pvalue_weight - pvalue[i].pvalue_weight_dist[j];
	   if(pvalue[i].pvalue_weight_dist[j] < epsilon) { pvalue[i].pvalue_weight_dist[j]=0;} 
	}
    }
  
 ////////Step 5:Computes coefficient  (products of weights) 
   
  Log_Coeff=0.0;
  for(i=0;i<int(pvalue.size());i++)  
    { 
      for(j=0;j<int(pvalue[i].pvalue_weight_dist.size()); j++)
	{
	  Log_Coeff = Log_Coeff + log(pvalue[i].pvalue_weight - pvalue[i].pvalue_weight_dist[j]);  
	}
    }
 ////////Step 7:Computes combine P-value 
 
  //cout("Log_Tau = %e\t Coeff = %e \n",Log_Tau,Log_Coeff);
  final_pvalue=Combine_Weighted_Pvalue_CDF(pvalue,Log_Tau,Log_Coeff);
  if(final_pvalue > 1.0) { final_pvalue = 1.0;} 
  
  return final_pvalue;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::Combine_Weighted_Pvalue_CDF(vector<struct CCW_COINED_PVALUE> pvalue,double Log_Tau,double Log_Coeff)
{
  int k,k1,j;
  int g[4]={0,2,3,4};
  double sum_pvalue=0.0;
  double f_pvalue,delta_k,delta_k1,approx;
  
  approx = Combine_Weighted_Pvalue_Approx(pvalue,Log_Tau);
  //cout("approx = %e\n",approx); 
  f_pvalue = exp( Log_Coeff + log(approx)); 
  if (f_pvalue == 0) { f_pvalue = DBL_MIN;}
  //cout("f_pvalue = %e \t approx = %e\n",f_pvalue,approx);  

  //Computes single summation terms correction 
  for(j=1;j<4;j++)
    {
      //first approximation
      sum_pvalue=0;
      for(k=0;k<int(pvalue.size());k++)
	{
	  //cout<<"k="<<k<<endl;
	  pvalue[k].pvalue_degeneracy = pvalue[k].pvalue_degeneracy +  g[j];
	  if(g[j] == 4)
            {
	      if( pvalue[k].pvalue_weight_dist.size() > 0)
	      {
		  delta_k = Combine_Weight_Pvalue_Delta_Weight(k,g[j],pvalue);
                  delta_k1 = Combine_Weight_Pvalue_Delta_Weight(k,2,pvalue); 
		   if(delta_k != 0 || delta_k1 != 0)
		    {
		      sum_pvalue= sum_pvalue + delta_k + 0.5*pow(delta_k1,2.0)*exp(Log_Coeff+log(Combine_Weighted_Pvalue_Approx(pvalue,Log_Tau)));
		      //cout<<"F1"<<k<<"="<<sum_pvalue<<endl;
		     }
		  }
	    }
           else
	     {
	       if( pvalue[k].pvalue_weight_dist.size() > 0)
	         {
		   delta_k = Combine_Weight_Pvalue_Delta_Weight(k,g[j],pvalue);
		    if( delta_k != 0)
		    {
		      sum_pvalue= sum_pvalue + delta_k*exp(Log_Coeff+log(Combine_Weighted_Pvalue_Approx(pvalue,Log_Tau)));
		      //cout<<"F2"<<k<<"="<<"delta_k="<<delta_k<<"\t"<<Combine_Weighted_Pvalue_Approx(pvalue,Log_Tau)<<endl;
		      }
		  }
	     }
	   pvalue[k].pvalue_degeneracy = pvalue[k].pvalue_degeneracy - g[j];
	}
          f_pvalue = f_pvalue + sum_pvalue;
	  // cout<<"Single summation P-value correction="<<f_pvalue<<endl;
     }
 
  //Computes double summation terms correction
    sum_pvalue=0;
     for(k=0;k<int(pvalue.size());k++)
	 {
            delta_k = Combine_Weight_Pvalue_Delta_Weight(k,g[1],pvalue); 
            if(delta_k != 0)
	      {
		for(k1=0;k1<int(pvalue.size());k1++)
		  {
		    if(k1 != k)
		      {
			pvalue[k].pvalue_degeneracy = pvalue[k].pvalue_degeneracy   +  g[1];
			pvalue[k1].pvalue_degeneracy = pvalue[k1].pvalue_degeneracy + g[1];
		 
			if((pvalue[k1].pvalue_weight_dist.size() > 0))
			  {
			    delta_k1 = Combine_Weight_Pvalue_Delta_Weight(k1,g[1],pvalue);
			    if(delta_k1 != 0)
			      {
				sum_pvalue= sum_pvalue + 0.5*delta_k1*delta_k*exp(Log_Coeff+log(Combine_Weighted_Pvalue_Approx(pvalue,Log_Tau)));
			      }
			  }
			pvalue[k].pvalue_degeneracy = pvalue[k].pvalue_degeneracy   -   g[1];
			pvalue[k1].pvalue_degeneracy = pvalue[k1].pvalue_degeneracy -   g[1];
		      }		
		  }
	      }  
	 }  
	 
     f_pvalue = f_pvalue + sum_pvalue;
     //cout<<"Double summation P-value correction="<<f_pvalue<<endl;    

  return f_pvalue;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Computes coefficient for correction terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double inline CoinedPValues::Combine_Weight_Pvalue_Delta_Weight(int k,int g,vector<struct CCW_COINED_PVALUE> pvalue)
{
  int j;
  double sum_delta=0.0; 

  for(j=0;j<int(pvalue[k].pvalue_weight_dist.size());j++)
    {
       
      sum_delta = sum_delta + pow(-1.0*pvalue[k].pvalue_weight_dist[j],g); 
      // cout("delta = k=%d g=%d j=%d   %e %e\n",k,g,j,sum_delta,pvalue[k].pvalue_weight_dist[j]);
    }

  return sum_delta/g;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Computing combined weighted P-value
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::Combine_Weighted_Pvalue_Approx(vector<struct CCW_COINED_PVALUE> pvalue,double Log_Tau)
{
  int i,k;
  int j,g,s,m;
  double product=1;
  double sum1,sum2,sum3,x;
  double epsilon = 1e-16;
  vector<int> index;
  vector<double> sum_positive;
  vector<double> sum_negative; 
  
  product = Log_Tau;
  m=pvalue.size()-1;
  sum3=sum2=sum1=0;
  for(k=0;k<=m;k++)
    {
      index.clear();  
      index.push_back(k);
      for(s=0;s<=m;s++) { if(s != k) { index.push_back(s);}}
	 
      sum2=0;   j=0; 
      if( m > 0)
	{
	  for(g=0;g<=pvalue[k].pvalue_degeneracy-1;g++) 
	    {
	      sum1=0; j=0;
	      for(s=0;s<=g;s++)
		{ 
		  x =  exp(-pvalue[k].pvalue_weight*product);
		  if(  x != 0)
		    {
		      sum1 = sum1 + x*pow(pvalue[k].pvalue_weight*product,s)/factorial[s];
		    }
		}

	      sum1=sum1/(pow(pvalue[k].pvalue_weight,g+1.0)); 
		
	      if(sum1 != 0)
		{
		  // sum2 = sum2 +  sum1*Recursion_Combine_Weighted_Pvalue_CDF(m,k,pvalue[k].pvalue_degeneracy -1-g,j,index,pvalue);
		  sum2 = sum1*Recursion_Combine_Weighted_Pvalue_CDF(m,k,pvalue[k].pvalue_degeneracy -1-g,j,index,pvalue);
		  if( sum2 < 0) { sum_negative.push_back(fabs(sum2));}
		  else { sum_positive.push_back(fabs(sum2));}
		}
	      //cout<<"1::s1="<<sum1<<"\t"<<"s2'="<<sum2<<endl;
	    }
	}
      else
	{
	  sum1=0;
	  g = pvalue[k].pvalue_degeneracy-1; 
	  for(s=0;s<=g;s++)
	    {
	      x =  exp(-pvalue[k].pvalue_weight*product);
	      if(  x != 0)
		{
		  sum1 = sum1 + x*pow(pvalue[k].pvalue_weight*product,s)/factorial[s];
		}
	    }
	     
	  // sum2 = sum2 + sum1;  
	  if( sum1 < 0) { sum_negative.push_back(fabs(sum1));}
	  else { sum_positive.push_back(fabs(sum1));}
	  //cout<<"2::s1="<<sum1<<"\t"<<"s2="<<sum2<<endl;
	}
    }

  sum1=sum2=0;
  for(i=0;i<int(sum_positive.size());i++) //Summing positive contribution 
    { 
      sum1 = sum1 + sum_positive[i];
      // cout<<"S-positive="<<sum_positive[i]<<"\t"<<sum1<<endl;
    }
  for(i=0;i<int(sum_negative.size());i++) //Summing negative contribution 
    { 
      sum2 = sum2 + sum_negative[i];  
      // cout<<"S-negative="<<sum_negative[i]<<"\t"<<sum2<<endl;  
    }
     
   if( sum1!=0)
   {
     // cout<<sum1<<"\t"<<sum2<<"\t"<<sum1-sum2<<endl;
       sum3 = sum1-sum2;
       if(sum3 <= 0)
 	{
 	  sum3=sum1*epsilon;
 	}
       // cout<<"Sum3 = "<<sum3<<endl;
     }
    else
     {
       sum3 = DBL_MIN;
     }

  return sum3;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::Recursion_Combine_Weighted_Pvalue_CDF(int m,int k,int n,int j,vector<int> index,vector<struct CCW_COINED_PVALUE> pvalue)
{
  int Gj; 
  double sum=0;
  double coef;

  j=j+1;

  if(j < m )
    { 
      for(Gj=0;Gj<=n;Gj++)
         {
	    coef =  factorial[pvalue[index[j]].pvalue_degeneracy +Gj-1]/(factorial[pvalue[index[j]].pvalue_degeneracy-1]*factorial[Gj])*pow(-1.0,Gj)/pow(pvalue[index[j]].pvalue_weight- pvalue[k].pvalue_weight,pvalue[index[j]].pvalue_degeneracy + Gj);

	   sum = sum + coef* Recursion_Combine_Weighted_Pvalue_CDF(m,k,n-Gj,j,index,pvalue);

	 }
    }
  else if(j==m)
    {
       Gj = n;
       coef = (factorial[ pvalue[index[j]].pvalue_degeneracy + Gj-1]/(factorial[pvalue[index[j]].pvalue_degeneracy - 1]*factorial[Gj]))*pow(-1.0,Gj)/pow(pvalue[index[j]].pvalue_weight - pvalue[k].pvalue_weight,pvalue[index[j]].pvalue_degeneracy + Gj);
    
       sum =sum + coef;
       
      
    }
  return sum;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*The function below computes the cummulative distribution function for the combined weighted P-value*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Mathai(vector<struct COINED_PVALUE> v_pvalue)
{
  double pvalue;
  pvalue = CoinedPvalues_Log_Mathai(v_pvalue);
  return exp(pvalue);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*The function below computes the cummulative distribution function for the combined weighted P-value*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_Mathai(vector<struct COINED_PVALUE> v_pvalue)
{
	printf("fuck yeah!\n");
  int k;
  int i,j,g,s,n,m;
  double product=0;
  double epsilon = 1e-6;
  double max_pvalue,log_p,log_pmax;
  double sum1,sum2,sum3,sr;
  double coef=1;
  vector<double> nd_pvalue;
  vector<double> nd_eigenvalue;
  vector<int> nd_degeneracy; 
  vector<int> index;
 
  Check_Pvalue_Range(v_pvalue); //Checking P-value between [0,1]

  //Removing P-value degeneracy
  for(i=0;i<int(v_pvalue.size());i++)
    {    
      if( v_pvalue[i].pvalue > 0 && v_pvalue[i].pvalue <=1 )
	{ 
	  nd_pvalue.push_back(log(v_pvalue[i].pvalue));
	  nd_eigenvalue.push_back(v_pvalue[i].pvalue_weight);     
	  nd_degeneracy.push_back(1); 
          n = nd_degeneracy.size()-1;
	  for(j=i+1;j<int(v_pvalue.size());j++)
	    {
	      
	      if( fabs(v_pvalue[i].pvalue_weight  -   v_pvalue[j].pvalue_weight) <= epsilon)
		{  
		  nd_pvalue[n] = nd_pvalue[n] + log(v_pvalue[j].pvalue);
		  nd_degeneracy[n] =   nd_degeneracy[n] + 1;
		  // cout<<"CoinedPValues::CoinedPvalues_Log_Mathai:P-value="<<v_pvalue[i].pvalue<<", weight="<<v_pvalue[i].pvalue_weight<<"\t"<<v_pvalue[j].pvalue<<", weight="<<v_pvalue[j].pvalue_weight<<", i="<<i<<", j="<<j<<", degeneracy="<<nd_degeneracy[n]<<endl;
		  v_pvalue[j].pvalue = -1;
		}  
	    }     
	}
    }

  //Computing P-value producy
  product = 0;
  //cout<<"CoinedPValues::CoinedPvalues_Log_Mathai:Number of pvalues="<<int(nd_eigenvalue.size())<<endl;
  for(i=0;i<int(nd_eigenvalue.size());i++)
    {
     
      product= product -  nd_eigenvalue[i]*nd_pvalue[i];
      nd_eigenvalue[i] = 1.0/nd_eigenvalue[i];
    }
 
  //Computing combine weight P-value CDF
     max_pvalue=-1e9;
     m=nd_pvalue.size()-1;
     sum3=sum2=sum1=0;
     for(k=0;k<=m;k++)
       {

	 n=nd_degeneracy[k]-1;
	 index.clear();  
	 index.push_back(k);
	 for(s=0;s<=m;s++) { if(s != k) { index.push_back(s);}}
         coef = coef*pow(nd_eigenvalue[k],nd_degeneracy[k]);
      
	 sum2=0;   j=0; 
	 if( m > 0)
	   {
	     // cout<<"CoinedPValues::CoinedPvalues_Log_Mathai:P-value degeneracy="<<nd_degeneracy[k]<<", weight="<<1.0/nd_eigenvalue[k]<<endl;
	     for(g=0;g<=nd_degeneracy[k]-1;g++) 
	       {
		 
		 sum1=0; j=0;  log_pmax=-1e9;
		 for(s=0;s<=g;s++)
		   {
                     log_p = -nd_eigenvalue[k]*product + s*log(nd_eigenvalue[k]*product) - FACTR_LN(s);
		     sum1 = sum1  + exp(log_p);
                     if( log_p > log_pmax) { log_pmax = log_p;}
		   }
                 
		 sum1= log(sum1) - (g+1)*log(nd_eigenvalue[k]);   
                 sr = Recursion_Combine_Weighted_Pvalue_CDF_Mathai(m,k,nd_degeneracy[k]-1-g,j,nd_eigenvalue,nd_degeneracy,index); 
		 sum2 = sum2 +  exp(sum1)*sr;

                 log_p = log_pmax; 
                 log_p = log_p - (g+1)*log(nd_eigenvalue[k]);
                 log_p= log_p + log(sr);
                 if(log_p > max_pvalue) { max_pvalue = log_p;}
	       }
	   }
	 else
	   {
	     sum1=0;
	     g = nd_degeneracy[k]-1; 
	     for(s=0;s<=g;s++)
	       { 
                 log_p = -nd_eigenvalue[k]*product + s*log(nd_eigenvalue[k]*product) - FACTR_LN(s);  
		 sum1 = sum1  + exp(log_p);
                 if( log_p > max_pvalue) { max_pvalue = log_p;}
	       }
	     sum2 = sum2 +  sum1;  
	   }

	 sum3=sum3+sum2;
	 //cout("Mathai P-value = %e\n",coef*sum2);
       }

   if( coef*sum3 > 0) { return log(coef*sum3);}
   else return  log(coef) + max_pvalue;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::Recursion_Combine_Weighted_Pvalue_CDF_Mathai(int m,int k,int n,int j,vector<double> nd_eigenvalue,vector<int> nd_degeneracy,vector<int> index)
{
  int Gj; 
  double sum=0;
  double coef,diff;

  j=j+1;

  if(j < m )
    { 
      for(Gj=0;Gj<=n;Gj++)
         {
           diff = nd_eigenvalue[index[j]]-nd_eigenvalue[k];   

           if(diff > 0)
	     { 
	       coef =  pow(-1.0,Gj)*exp(FACTR_LN(nd_degeneracy[index[j]]+Gj-1) - FACTR_LN(nd_degeneracy[index[j]]-1) - FACTR_LN(Gj) -  (nd_degeneracy[index[j]]+ Gj)*log(diff));
	     }
           else
	     {
	       coef =  pow(-1.0,Gj-(nd_degeneracy[index[j]]+ Gj))*exp(FACTR_LN(nd_degeneracy[index[j]]+Gj-1) - FACTR_LN(nd_degeneracy[index[j]]-1) - FACTR_LN(Gj) -  (nd_degeneracy[index[j]]+ Gj)*log(fabs(diff)));
	     }

	   sum = sum + coef*Recursion_Combine_Weighted_Pvalue_CDF_Mathai(m,k,n-Gj,j,nd_eigenvalue,nd_degeneracy,index);
	 }
    }
  else if(j==m)
    {
       Gj = n;
      
        diff = nd_eigenvalue[index[j]]-nd_eigenvalue[k];   

           if(diff > 0)
	     { 
	       coef =  pow(-1.0,Gj)*exp(FACTR_LN(nd_degeneracy[index[j]]+Gj-1) - FACTR_LN(nd_degeneracy[index[j]]-1) - FACTR_LN(Gj) -  (nd_degeneracy[index[j]]+ Gj)*log(diff));
	     }
           else
	     {
	       coef =  pow(-1.0,Gj-(nd_degeneracy[index[j]]+ Gj))*exp(FACTR_LN(nd_degeneracy[index[j]]+Gj-1) - FACTR_LN(nd_degeneracy[index[j]]-1) - FACTR_LN(Gj) -  (nd_degeneracy[index[j]]+ Gj)*log(fabs(diff)));
	     }
       sum =sum + coef;
    }
  return sum;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combining Truncate-Weighted-P-value product using Good's Formula.*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Good( vector<struct COINED_PVALUE> v_pvalue)
{
  double pvalue; 
  pvalue = CoinedPvalues_Log_Good(v_pvalue);
  return exp(pvalue);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combining Truncate-Weighted-P-value product using Good's Formula.*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_Good( vector<struct COINED_PVALUE> pvalue)
{
  int i,j,n,log_sum_imag;
  long double p, product,log_sum_weights;
  long double log_p,p_tau;
  double tau = TAU;
  
  //Checking P-value between [0,1]
  Check_Pvalue_Range(pvalue); 

  //Combining P-value with identical weights
   Remove_Degenerate_Weights(pvalue);

  n=pvalue.size();
  if(n > 0)
    {
      //computing product of P-value in logarithmic scale
      product=0;
      for(i=0;i<n;i++) 
	{  
	  product =  product + pvalue[i].pvalue_weight*log(pvalue[i].pvalue);
	}
  
      p=0; 
      for(i=0;i<n;i++) 
	{
	  if(tau < 1)//Used only when combining truncate P-values 
	    {
	      p_tau=0;
	      for(j=0;j<n;j++) 
		{
		  p_tau = p_tau  + (pvalue[i].pvalue_weight - pvalue[j].pvalue_weight);
		}
	      //cout<<i<<"\t"<<p_tau<<"\t"<<tau<<"\t"<<n<<"\t"<<p_tau - n<<"\t"<<p_tau/pvalue[i].pvalue_weight - n<<endl;
	      p_tau = p_tau/pvalue[i].pvalue_weight - n;  
	      p_tau = p_tau*log(tau); 
	    }
	  else
	    {
	      p_tau=0;
	    }

	  if( n > 1)
	    {
	      Good_Denominator(i,pvalue,log_sum_weights,log_sum_imag);
	      if( fmod(log_sum_imag,2.0) == 0)
		{
		  log_p =  product/pvalue[i].pvalue_weight +  (n-1)*log(pvalue[i].pvalue_weight) - log_sum_weights + p_tau; 
		  // cout<<exp(log_p)<<"\t"<<product<<"\t"<<pvalue[i].pvalue_weight<<"\t"<<(n-1)*log(pvalue[i].pvalue_weight)<<"\t"<<log_sum_weights<<"\t"<<p_tau<<endl;
		  p = p + exp(log_p);
		 
		}
	      else
		{
		  log_p =  product/pvalue[i].pvalue_weight +  (n-1)*log(pvalue[i].pvalue_weight) - log_sum_weights + p_tau;
		  // cout<<exp(log_p)<<"\t"<<product<<"\t"<<pvalue[i].pvalue_weight<<"\t"<<(n-1)*log(pvalue[i].pvalue_weight)<<"\t"<<log_sum_weights<<"\t"<<p_tau<<endl;
		  p = p - exp(log_p); 
		}
	    }
	  else
	    {
	      log_p =  product  + p_tau; 
	      //cout<<exp(log_p)<<"\t"<<product<<"\t"<<pvalue[i].pvalue_weight<<"\t"<<(n-1)*log(pvalue[i].pvalue_weight)<<"\t"<<log_sum_weights<<"\t"<<p_tau<<endl;
	      p = p + exp(log_p);  
	    }
	}

      if(p <=1) 
	{
	  if( p > 0 ) return log(p);
	  else return log(DBL_MIN);
	 }   
       else return 0;
    }
  else return 0;
}
void CoinedPValues::Good_Denominator(int v, vector<struct COINED_PVALUE> &pvalue,long double &log_sum_weights, int &log_sum_imag)
{
  int i;
  long double diff;
  //computing log of weight
  log_sum_weights=log_sum_imag=0;
  for(i=0;i<int( pvalue.size());i++)
    {
      if( i != v)
	{
	  diff = pvalue[v].pvalue_weight - pvalue[i].pvalue_weight;
	  log_sum_weights = log_sum_weights + log(fabs(diff));
          if(diff < 0) { log_sum_imag = log_sum_imag +1;}
	}
    }  //cout<<v<<"Good demoniator log sum = "<<log_sum_weights<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combines weighted P-value using Fisher's Formula*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Fisher( vector<struct COINED_PVALUE> v_pvalue)
{
  double pvalue;
  pvalue = CoinedPvalues_Log_Fisher(v_pvalue);
  return exp(pvalue);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combines weighted P-value using Fisher's Formula*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_Fisher( vector<struct COINED_PVALUE> pvalue)
{
  int i,n;
  double product,p;
  double log_p,p_tau;
  Check_Pvalue_Range(pvalue); //Checking P-value between [0,1]

  n=pvalue.size();
  product=0;
  for(i=0;i<n;i++)  
    { 
      product =  product + log(pvalue[i].pvalue); 
    }

  p=0; 
  if(product < 0)
    {
      if(TAU < 1) { p_tau = n*log(TAU);}
      else{ p_tau = 0;}
      for(i=0;i<n;i++) 
	{
	  log_p = (product - p_tau + i*log(fabs(product) + p_tau) - FACTR_LN(i));
	  p = p + exp(log_p);
	}
      if( p > 0) { return log(p);}
      else return log(DBL_MIN); 
    }
  else return 0;
}
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combines P-value using Fisher's Formula with Tau==1*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Fisher_Tau1(vector<struct COINED_PVALUE> pvalue)
{
  int i,n;
  double product,p,tau;
  double log_p,p_tau;
  Check_Pvalue_Range(pvalue); //Checking P-value between [0,1]
 
  n=pvalue.size();
  product=0;
  tau=1;
  for(i=0;i<n;i++)  
    { 
      product =  product + log(pvalue[i].pvalue); 
       if( pvalue[i].pvalue > 1) 
	 { 
	   cerr<<"CoinedPValues::Error:CoinedPvalues_Log_Fisher:( pvalue[i].pvalue > 1)"<<endl; exit(1);
	 }    
    }
  if(tau < 1) { p_tau = n*log(tau);}
  else{ p_tau = 0;}

  p=0;
  if( n > 0)
    {
      for(i=0;i<n;i++) 
	{
	  log_p = (product - p_tau + i*log(fabs(product) + p_tau) - FACTR_LN(i));
	  p = p + exp(log_p);
	}
      if(p > 0) return p;
      else return DBL_MIN;
    }
  else return  1;
}
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combines P-value using Fisher's Formula with Tau==1*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Fisher_Log_Tau1(vector<struct COINED_PVALUE> pvalue)
{
  int i,n;
  double product,p,tau;
  double log_p,p_tau;
  Check_Pvalue_Range(pvalue); //Checking P-value between [0,1]
 
  Remove_Degenerate_Weights(pvalue); //Combining degenerate P-values

  n=pvalue.size();
  product=0;
  tau=1;
  for(i=0;i<n;i++)  
    { 
      product =  product + log(pvalue[i].pvalue); 
       if( pvalue[i].pvalue > 1) 
	 { 
	   cerr<<"CoinedPValues::Error:CoinedPvalues_Log_Fisher:( pvalue[i].pvalue > TAU)"<<endl; exit(1);
	 }    
    }
  if(tau < 1) { p_tau = n*log(tau);}
  else{ p_tau = 0;}

  p=0; 
  if( n > 0)
    {
      for(i=0;i<n;i++) 
	{
	  if(fabs(product) + p_tau <= 0) { cerr<<"CoinedPValues::CoinedPvalues_Log_Fisher: fabs(product) + p_tau <= 0)"<<endl; exit(1);}
	  log_p = (product - p_tau + i*log(fabs(product) + p_tau) - FACTR_LN(i));
	  p = p + exp(log_p);
	}
      if( p > 0) { return log(p);}
      else return log(DBL_MIN);
    }
  else return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Return the combine P-value using Bhoj approximation*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Bhoj(vector<struct COINED_PVALUE> v_pvalue)
{
  double pvalue; 
  pvalue = CoinedPvalues_Log_Bhoj(v_pvalue);
  return exp(pvalue);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Return the combine P-value using Bhoj approximation*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_Bhoj(vector<struct COINED_PVALUE> pvalue)
{
  int i;
  double product,normalization_w,p,weight;
  double log_p;

  Check_Pvalue_Range(pvalue); //Checking P-value between [0,1]
  product=normalization_w=0;
  
  for(i=0;i<int( pvalue.size());i++)  
    { 
      normalization_w = normalization_w + pvalue[i].pvalue_weight;
      product =  product +  pvalue[i].pvalue_weight*log(pvalue[i].pvalue); 
    }
  product=-2*product/normalization_w;

  p=0;
  if( pvalue.size() > 0)
    {
      for(i=0;i<int(pvalue.size());i++)  
	{
	  weight = (pvalue[i].pvalue_weight/normalization_w);
	  log_p = ln_gammq(1.0/weight,product/(2*weight)); //log of Q(a,x)
	  log_p = log(weight) + log_p; 
	  p = p + exp(log_p);
	}
      if( p > 0) { return log(p);}
      else return log(DBL_MIN);
    }
  else return  0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combined Weighted P-values (WP) using gamma distribution*/ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_WP_RTGD(vector<struct COINED_PVALUE> v_pvalue)
{
  double pvalue; 
  pvalue = CoinedPvalues_Log_WP_RTGD(v_pvalue);
  return exp(pvalue); 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combined Weighted P-values (WP) ratio of two gamma distribution RTGD*/ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_WP_RTGD(vector<struct COINED_PVALUE> pvalue)
{
  int i;
  double normalize_weight;
  double weight,B1,B2,B3;
  double alpha_1,alpha_2,gamma;
  double X,p;

  Check_Pvalue_Range(pvalue); //Checking P-value between [0,1]

  //Normalizing weights
  normalize_weight = 0;
  for(i=0;i<int(pvalue.size());i++) { normalize_weight = normalize_weight + pvalue[i].pvalue_weight; }
 
  B1=B2=B3=0; 
  X = 0;
  for(i=0;i<int(pvalue.size());i++)
    {
      weight = pvalue[i].pvalue_weight/normalize_weight; 
      X = X - weight*log(pvalue[i].pvalue);
      B1 = B1 + weight;
      B2 = B2 + pow(weight,2.0);
      B3 = B3 + pow(weight,3.0);
    }
  //Computing distribution moments
  alpha_1 = B1*(B1*B1*B2 + 2*B1*B3 - B2*B2)/(2*B1*B2*B2 - B1*B1*B3 + B2*B3);
  alpha_2 = (B1*B1*B2 + 3*B1*B3 - 2*B2*B2)/(B1*B3-B2*B2);
  gamma = (2*B1*B2*B2 - B1*B1*B3 + B2*B3)/(B1*B3 - B2*B2);
  
  //cout<<"B1 = "<<B1<<endl;
  //cout<<"B2 = "<<B2<<endl;
  //cout<<"B3 = "<<B3<<endl;
  //cout<< alpha_1<<"\t"<< alpha_2<<"\t"<<gamma<<endl;
  
  if(  alpha_1 > 0 &&  alpha_2 > 0)
    {
      //Applying distribution transformation
      X = (X/gamma)/(X/gamma + 1);   
      p = Pvalue_Beta_Pdf(alpha_1, alpha_2,X);
    }
  else
    {
      alpha_1 = B1*B1/B2;   
      alpha_2 = B1/B2;
      p = ln_gammq(alpha_1,alpha_2*X);  
    }
  return p;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combined Weighted P-values (WP) using gamma distribution*/ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_WP_GD(vector<struct COINED_PVALUE> v_pvalue)
{
  double pvalue; 
  pvalue = CoinedPvalues_Log_WP_GD(v_pvalue);
  return exp(pvalue); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Combined Weighted P-values (WP) using gamma distribution*/ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_WP_GD(vector<struct COINED_PVALUE> pvalue)
{
  int i;
  double normalize_weight;
  double weight,B1,B2;
  double alpha_1,alpha_2;
  double X,p;

  Check_Pvalue_Range(pvalue); //Checking P-value between [0,1]

  //Normalizing weights
  normalize_weight = 0;
  for(i=0;i<int(pvalue.size());i++) { normalize_weight = normalize_weight + pvalue[i].pvalue_weight; }

  B1=B2=0; 
  X = 0;
  for(i=0;i<int(pvalue.size());i++)
    {
      weight = pvalue[i].pvalue_weight/normalize_weight; 
      X = X - weight*log(pvalue[i].pvalue);
      B1 = B1 + weight;
      B2 = B2 + pow(weight,2.0);
    }
  //Computing distribution moments 
  alpha_1 = B1*B1/B2;   
  alpha_2 = B1/B2;
  p = ln_gammq(alpha_1,alpha_2*X); 
  return p;
}
/////////////////////////////////////////////////////////////////////////////////////
/*Returns the weighted combine truncated P-value using the effective dimension equal to the sum of weights*/
////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_TPM_EFFD_Fisher(vector<struct COINED_PVALUE> pvalue)
{
  int i,s,L;
  double t,c_pvalue,w,log_w;
  double tau;

  L = int(pvalue.size());
  log_w=w=0; 
  for(i=0;i<L;i++) 
    { 
      log_w  = log_w + pvalue[i].pvalue_weight*log(pvalue[i].pvalue);  //product of weighted p-values
      w = w + pvalue[i].pvalue_weight; //effective dimension
    } 

  tau = w*log(TAU);  
  if( w < 1) { L = 1;}
  else { L = int(w+0.5);}
  w=log_w;
  c_pvalue=0; 
 
  if(w < tau)
   {
     t=0;
     for(s=0;s<=L-1;s++)
       {
	 t = t +  exp(s*log(tau - w)-FACTR_LN(s));
       }
     if( t > DBL_MAX) { t = DBL_MAX;}
     if( t > 0)
       {
	 c_pvalue = w - tau + log(t);
       }
     else 
       {
	c_pvalue=0; 
       }
   }
 else
   {
     c_pvalue=0;
   }  
 return c_pvalue;
}
/////////////////////////////////////////////////////////////////////////////////////
/*Returns the weighted combine truncated P-value using the effective dimension equal to the sum of weights*/
////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_Log_TPM_EFFD_Fisher_Tau1(vector<struct COINED_PVALUE> pvalue)
{
  int i,s,L;
  double t,c_pvalue,w,log_w;
  double tau;

  L = int(pvalue.size());
  log_w=w=0; 
  for(i=0;i<L;i++) 
    { 
      log_w  = log_w + pvalue[i].pvalue_weight*log(pvalue[i].pvalue);  //product of weighted p-values
      w = w + pvalue[i].pvalue_weight; //effective dimension
    } 

  tau = 0;  // tau=w*log(TAU); Tau=1
  if( w < 1) { L=1;}
  else  L=int(w+0.5);
  w=log_w;
  c_pvalue=0; 
 
  if(w < tau)
    {
      t=0;
      for(s=0;s<=L-1;s++)
	{
	  t = t +  exp(s*log(tau - w)-FACTR_LN(s));
	}
      if( t > DBL_MAX) { t = DBL_MAX;}
      if( t > 0)
	{
	  c_pvalue = w - tau + log(t);
	}
      else 
	{
	  c_pvalue=0; 
	}
    }
  else
    {
      c_pvalue=0;
    }  
  return c_pvalue;
}
/////////////////////////////////////////////////////////////////////////////////////
/*Returns the unweighted combine truncated P-value*/
//Parameter:( pvalue_vector)
////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::CoinedPvalues_TPM_Zaykin(vector<struct COINED_PVALUE> pvalue)
{
  int L,i;
  int k,s;
  double b,t,w,c_pvalue,log_w;

  L = int(pvalue.size());
  log_w=0;
  for(i=0;i<L;i++) { log_w  = log_w + log(pvalue[i].pvalue); } //product of p-values

  //cout<<"CoinedPValues::CoinedPvalues_TPM_Zaykin: log(P-value) = "<<log_w<<endl;

  w = log_w;
  if( L >= 1)
    {
      c_pvalue=0;
      for(k=1;k<=L;k++)
	{
	  b =  FACTR_LN(L) - FACTR_LN(k) - FACTR_LN(L-k);      
	  if(w <= k*log(TAU))
	    {
	      t=0;
	      for(s=0;s<=k-1;s++) { t = t +  exp(s*log(k*log(TAU) - w)-FACTR_LN(s));}
         
	      if( t > 0){ t= w+log(t)+ b +(L-k)*log(1-TAU);}
	    }
	  else
	    {
	      t = k*log(TAU);
	      t= t + b + (L-k)*log(1-TAU);	     
	    }	   
	  if( t < 0)
	    {
	      c_pvalue = c_pvalue + exp(t);
	      //condition to terminate loop;
	      if(c_pvalue > 0) { if( fabs(t - log(c_pvalue)) > 10 ) { k=L;} }
	    } 
	}
    }
  else
    {
      c_pvalue=1;
    }

  //cout<<"CoinedPValues::CoinedPvalues_TPM_Zaykin: P-value = "<<c_pvalue<<endl;

  return c_pvalue;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Checking user input P-value is between [0,1]*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoinedPValues::Check_Pvalue_Range(vector<struct COINED_PVALUE> &v_pvalue)
{
  int i;

  for(i=0;i<int(v_pvalue.size());i++)  
    {  
      if( v_pvalue[i].pvalue == 0)
	{
	  v_pvalue[i].pvalue = DBL_MIN;
	}
      else if( v_pvalue[i].pvalue < 0 || v_pvalue[i].pvalue > 1)
	{
          cerr<<"P-values entered is not between [0,1] "<<v_pvalue[i].pvalue<<"\t"<<endl;
	  exit(1);
	}
    }
}
///////////////////////////////////////////////////////////////////
/*Returns  the incomplete cumulative beta pdf function integral.
Returns the incomplete beta function Ix(a; b).*/
///////////////////////////////////////////////////////////////////
double CoinedPValues::Pvalue_Beta_Pdf(double a, double b, double x)
{
double bt;

 if (x < 0.0 || x > 1.0) cerr<<"Bad x in routine betai"<<endl; 
 if (x == 0.0)  bt=0.0;
 else if (x==1) bt=1.0; 
 else  // Factors in front of the continued fraction.
   bt=GAMMA_LN(a+b)-GAMMA_LN(a)-GAMMA_LN(b)+a*log(x)+b*log(1.0-x);
 if (x < (a+1.0)/(a+b+2.0)) //Use continued fraction directly.
   {
     if( (bt + log(betacf(a,b,x)/a)) < 1)
       {
	 return log(1.0 - (bt + log(betacf(a,b,x)/a)));
       }
     else
       {
	 return -1e9;   
       }
   }
 else //Use continued fraction after making the symreturn
   return  bt + log(betacf(b,a,1.0-x)/b); //metry transformation.
}
////////////////////////////////////////////////////////////////////
/*Used by betai: Evaluates continued fraction for incomplete beta 
function by modied Lentz's  method (x5.2).*/
////////////////////////////////////////////////////////////////////
double CoinedPValues::betacf(double a, double b, double x)
{
  int m,m2;
  int ITMAX = 100;
  double aa,c,d,del,h,qab,qam,qap;
  double FPMIN = DBL_MIN*10;    //Number near the smallest representable
  double  EPS = 1e-7;     //Relative accuracy  
  qab=a+b; // These q's will be used in factors that occur
  qap=a+1.0; //in the coecients (6.4.6).
  qam=a-1.0;
  c=1.0;// First step of Lentz's method.
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=ITMAX;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d; //One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; //Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d=FPMIN; 
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break; // Are we done?
  }
  if (m > ITMAX) cerr<<"a or b too big, or ITMAX too small in betacf"<<endl;
return h;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Returns ln(n!) */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::FACTR_LN(int n)
{
  static int ntop=1;
  static double a[5001]={0,0}; //Fill in table only as required.
  int j;
  if (n < 0) cerr<<"Negative factorial in routine factr"<<endl; 
  if (n > 5000) return GAMMA_LN(n+1.0);
  //Larger value than size of table is required. Actually, this big a value is going to overflow on many computers, but no harm in trying.
  while (ntop<=n) 
    { 
      j=ntop++;
      a[ntop]=a[j] + log(ntop);
    }
  return a[n];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.  
  GAMMA(n)= (n-1)!. Returns ln[gamma(xx)]*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::GAMMA_LN(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
   int j;
   y=x=xx;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.000000000190015;
   for (j=0;j<=5;j++) ser += cof[j]/++y; 
return -tmp+log(2.5066282746310005*ser/x);    //Returns the value ln[ÃƒÂƒ(xx)] for xx > 0.
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Return Q(a,x) the value of the complement of P(a,x)
that is Q(a,x) = 1 - P(a,x)*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CoinedPValues::ln_gammq(double a,double x)
{
  double gamser,gammcf,gln;

  if(x < 0.0 || a <= 0.0) { cout<<"Invalid arguments in gammq. a = "<<a<<",x="<<x<<endl; }
  if( x < (a+1)) 
    {
      ln_gser(&gamser,a,x,&gln);
      if( (1 - exp(gamser) == 0)) 
	return -1000000;
      else
	return log (1 - exp(gamser));
    }
  else
    {
      ln_gcf(&gammcf,a,x,&gln);
      return gammcf;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation
as gammcf. Also returns ln ??(a) as gln.*/
////////////////////////////////////////////////////////////////////
void CoinedPValues::ln_gcf(double *gammcf, double a, double x, double *gln)
{
  int i;
  int ITMAX=100;
  double an,b,c,d,del,h;
  double FPMIN = DBL_MIN*10;    //Number near the smallest representable
  double  EPS = DBL_EPSILON;     //Relative accuracy  

  *gln= GAMMA_LN(a);
  b=x+1.0-a; //Set up for evaluating continued fraction by modied Lentz's method (x5.2) with b0 = 0.
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;

  // cout<<"c = "<<c<<"\t"<<"d = "<<d<<"\t"<<"h = "<<"\t"<<h<<"\t"<<endl;
  // cout<<"FPMIN = "<<FPMIN<<"\t"<<"EPS = "<<EPS<<endl<<endl;
 
  for (i=1;i<=ITMAX;i++) 
    { 
      an = -i*(i-a);
      b += 2.0;
      d=an*d+b;
      if (fabs(d) < FPMIN) d=FPMIN;
      c=b+an/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;

      // cout<<"del ="<<del<<"\t"<<(-x+a*log(x)-(*gln))*h<<"\t"<<exp(-x+a*log(x)-(*gln))*h<<endl;

      if (fabs(del-1.0) < EPS) break;
    }
  if (i > ITMAX) cerr<<"a too large, ITMAX too small in gcf"<<endl;
  //*gammcf=exp(-x+a*log(x)-(*gln))*h; // Put factors in front.
     *gammcf=(-x+a*log(x)-(*gln))+log(h);  //log
   
  // cerr<<*gammcf<<endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////
/*Returns the incomplete gamma function P(a; x) evaluated by its series representation as gamser.
  Also returns ln ??(a) as gln.*/
////////////////////////////////////////////////////////////////////////////////////////////////
void CoinedPValues::ln_gser(double *gamser, double a, double x, double *gln)
{
   int n;
   int ITMAX=100;
   double sum,del,ap;
   double  EPS = DBL_EPSILON;     //Relative accuracy   

   *gln=GAMMA_LN(a);
   if (x <= 0.0) {
   if (x < 0.0) cerr<<"x less than 0 in routine gser"<<endl;
   *gamser=0.0;
   return;
   } else {
   ap=a;
   del=sum=1.0/a;
   for (n=1;n<=ITMAX;n++) {
   ++ap;
   del *= x/ap;
   sum += del;
   if (fabs(del) < fabs(sum)*EPS)
     {
       // *gamser=sum*exp(-x+a*log(x)-(*gln));
       *gamser = log(sum) + (-x+a*log(x)-(*gln));
       return;
     }
   }
   cerr<<"a too large, ITMAX too small in routine gser"<<endl;
return;
   }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Compute integer factorial terms
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoinedPValues::Factorial()
{
 int i;
 factorial[0]=1;
 for(i=1;i<=150;i++) { factorial[i]=i*factorial[i-1];}
 for(i=151;i<=5000;i++) { factorial[i]=sqrt(2*i*acos(-1.0))*pow(i/exp(1.0),1.0*i);}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
