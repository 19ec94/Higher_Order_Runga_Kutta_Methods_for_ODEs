
#include <iostream>
#include <cmath>
#include "functions.h"
#include <chrono>
#include <ctime>
#include <fstream>

using namespace std;

template<int order>
inline void RTOL(double *rtol_p);

int main()
{
 int *idid=NULL, *iwork=NULL, *spcrad=NULL,index;
 const int ns=250, neqn=ns*ns*2 ,nssq=ns *ns; 
 double *y=NULL,error=0.0;
 double *truey=NULL,ans,xx,yy;
 double t,tend,r;
 double h;
 double rtol,Atol;
 double *work=NULL;
 
 idid = new int;
 iwork = new int[10];
 spcrad = new int[2];
 y = new double[neqn];
 truey = new double[neqn];
 work = new double[1+5*neqn]; 
 FILE *sol, *sol_xy;
 
 sol = fopen("sol_serk5.txt","w");  //NEEDS TO BE CHANGED according to ORDER
 sol_xy =fopen("sol_serk5_xy.txt","w");


 t=0.0;
 tend=1.0; 
 ans=ns;

 for(int j=1; j<=ns; j++) 
 {
  yy=(j-1)/ans; 
  for(int i=1; i<=ns; i++)
  {
    y[ ((j-1)*ns+i )*2-2 ]=22.0*yy* pow ((1.0-yy),1.5);
  }
 }

 for(int i=1; i<=ns; i++) 
 {
  xx=(i-1)/ans; 
  for(int j=1; j<=ns; j++)
  {
    y[ ((j-1)*ns+i )*2-1 ]=27.0*xx* pow ((1.0-xx),1.5);
  }
 }

 RTOL<ORDER>(&rtol);
 Atol=rtol;  
 h=rtol;
 spcrad[0]=1;
 spcrad[1]=1;

 *idid=0;
 /********************************
  Initialize the Integration
**********************************/
 auto start=std::chrono::system_clock::now();
 ESERK(neqn, t, tend, h,   y, Atol,  spcrad, iwork,   work, idid);
 auto end=std::chrono::system_clock::now();
 std::chrono::duration<double> elapsed_seconds =end-start;
 std::time_t end_time = std::chrono::system_clock::to_time_t(end);
 cout<<"Time needed  "<<elapsed_seconds.count()<<endl;
 
 /***********************************
  * Print solution to files
  * *********************************/
 for(int j=1; j<=ns; j++) 
  {
  for(int i=1; i<=ns; i++)
  {
   fprintf(sol,"%d \t %.16f\n",i,y[ ((j-1)*ns+i )*2-2 ] );
  }
 }
  printf("\n");
  for(int i=1; i<=ns; i++) 
  {
  for(int j=1; j<=ns; j++)
  {
   fprintf(sol,"%d \t %.16f\n",j,y[ ((j-1)*ns+i )*2-1 ] );
  }
 }
 
 
 fclose(sol);
  
/*------------------------------------------------------
  Done.  Report some statistics.
------------------------------------------------------*/
/* ----- print statistics -----*/
     cout<< "Solution is tabulated in file sol_serk_xy.out"<<endl;
      cout<<"The value of IDID is "<<*idid<<endl;
      //cout<<"Max estimation of the spectral radius= "<<iwork[11]<<endl;
      //cout<<"Min estimation of the spectral radius= "<<iwork[12]<<endl;
      cout<<"Max number of stages used=" <<iwork[10]<<endl;
      cout<<"Number of f eval. for the spectr. radius="<<iwork[9]<<endl;
      cout<<"Number of f evaluations= "<<iwork[5]<<" steps= "<<iwork[7]+iwork[8]
            <<" accpt="<<iwork[7]<<" rejct= "<<iwork[8]<<endl;

/*------------------------------------------------------
  Done.  Compute the error and report some statistics.
------------------------------------------------------*/
/*   
   std::ifstream infile("sol_exact.out");
   int i=0;
   while (infile >> truey[i]){}
   infile.close();
   k=0;
    for (int j=1; j<=ns; j++)
    {
      for(int i=1 i<=ns; i++)
      {
        error = max(error,fabs(y[ ((j-1)*ns+i)*2-2] - truey[k] ));
        error = max(error,fabs(y[ ((j-1)*ns+i)*2-1] - truey [ns*ns+k] ) );
        k=k+1;
      }
    }
*/
  //delete  idid;
  //delete [] iwork;
  //delete [] spcrad;
  delete [] y;
  //delete [] truey;
 // delete [] work;
return 0;
}

template<int order>
inline void RTOL(double *rtol_p)
 {
   if(order ==4) {*rtol_p = pow(10.0,-3);}
   else          {*rtol_p = pow(10.0,-3);}
   return ;
 }


inline double rho(int neqn,double t, double *y)
{
 double alpha=0.25, nssq =neqn/2;
 return (8.0*nssq*alpha+2.0);
}



















