#include <iostream>
#include <cmath>
#include <algorithm>
#include "Constants1.h"
//SERKrho(neqn,t,f,yn,fn,iwork,hmax,work,sprad,idid)

using namespace std;

/*
int main()
{
 int neqn=25;
 int t =0; 
 long double *yn = new long double[neqn]; 
 long double *fn = new long double[neqn];
 long double *iwork = new long double[10];
 long double hmax = abs(1.48-t);
 long double *work = new long double[1+5*neqn];
 long double *sprad= new long double;
 int *idid= new int; 
*/
long double SERKrho(int neqn,long double t,long double *yn,long double *fn,int *iwork,long double hmax,long double *work,long double *sprad,int *idid)
{
 const int itmax=100;
 int i,iter,index,ptr5;
 long double uround=2.22L*pow(10,-16),sqrtu,ynrm;
 long double sigma=0,sigmal=0,dynrm=0,dfnrm=0,vnrm,small=0;
 long double *v=NULL,*fv=NULL;
 int nsteps=0,nfesig=0;
 
 
 sqrtu =sqrt(uround);
 nfesig = iwork[9];
 small = (long double)1.0/hmax;
 ptr5 = 4*neqn;
 v = new long double[neqn]; 
 fv = new long double[neqn];

 nsteps =iwork[7]+iwork[8];
 
 if(nsteps ==0)
  for(int i=0; i<neqn; i++){v[i]=fn[i];}
 else
  for(int i=0; i<neqn; i++){v[i]=work[ptr5+i];}
 
 ynrm=0.0L;
 vnrm=0.0L;
 
 for(int i=0; i<neqn; i++)
 {
   ynrm = ynrm +pow(yn[i],2);
   vnrm = vnrm +pow(v[i],2);
 }
 ynrm = sqrt(ynrm);
 vnrm = sqrt(vnrm);
 //cout<<"ynrm "<<ynrm<<endl;
 //cout<<"vnrm "<<vnrm<<endl;

 if((ynrm != 0.0) && (vnrm != 0.0))
 {
   dynrm = ynrm*sqrtu;
   for(int i=0; i<neqn; i++){v[i]=yn[i]+v[i]*(dynrm/vnrm);}
 }
 else if(ynrm !=0.0)
 {
   dynrm = ynrm*sqrtu;
   for(int i=0; i<neqn; i++){v[i]=yn[i]+yn[i]*sqrtu;}
 }
 else if(vnrm !=0.0)
 {
   dynrm = uround;
   for(int i=0; i<neqn; i++){v[i] = v[i]*(dynrm/vnrm);}
 }
 else
 {
   dynrm = uround;
   for(int i=0; i<neqn; i++){v[i] = dynrm;}
 }
 //cout<<"dynrm "<<dynrm<<endl;
 //Now iterate with a nonlinear power method
 sigma =0.0L;
 for(int iter=1; iter<=itmax; iter++)
 { //cout<<"Iter :"<<iter<<endl;
   fc(neqn,t,v,fv);
 //  for(int i=0; i<neqn; i++)cout<<fv[i]<<'\t';
  // cout<<'\n';
   nfesig = nfesig+1;
   dfnrm =0.0L;
   for(int i=0; i<neqn; i++){dfnrm = dfnrm + pow((fv[i]-fn[i]),2);}
   dfnrm = sqrt(dfnrm);
   //cout<<"dfnrm "<<dfnrm<<endl;
   sigmal = sigma;
   sigma = dfnrm/dynrm;
   *sprad = (12.0/10.0)*sigma;
   if((iter >= 2) && ( abs(sigma-sigmal) <= max(sigma,small)*0.01 ) )
   {
     for(int i=0; i<neqn; i++)
     {
       work[i]=yn[i];
       work[neqn+i]=fn[i];
       work[2*neqn+i]=v[i];
       work[3*neqn+i]=fv[i];
       work[ptr5+i]=v[i]-yn[i];       
     }
     iwork[9]=nfesig;
     //cout<<"nfesig "<<nfesig<<endl;
     //cout<<"sprad "<<*sprad<<endl;
     return 0;
   }
   if(dfnrm != 0.0)
   { 
     for(int i=0; i<neqn; i++)
     {
       v[i]=yn[i]+(fv[i]-fn[i])*(dynrm/dfnrm);
     }
   }
   else
   {
     index = 1+(iter%neqn);        //this would never hapen but asjjḱ , whether it should be 0 or 1
     v[index]=yn[index]-(v[index]-yn[index]);
   } 
  }
  iwork[9] = nfesig;
  *idid=6;


  delete [] v;
  delete [] fv;
  
 return  0;
}

