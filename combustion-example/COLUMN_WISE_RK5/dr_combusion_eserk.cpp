
#include <iostream>
#include <cmath>
#include <chrono>
#include <ctime>
#include "Constants1.h"
using namespace std;


int main()
{
 int i,j,*idid=NULL, *iwork=NULL, *spcrad=NULL;
 const int ns=100, neqn=ns*ns; 
 long double *y=NULL,error,*truey=NULL,ans,xx,yy;
 long double t,tend,r,h,rtol,Atol,*work=NULL;
 //External function has to be declared:

 long double etime, elapsed[2]; 
 int now[3];
 idid = new int;
 iwork = new int[10];
 spcrad = new int[2];
 y = new long double[neqn];
 truey = new long double[neqn];
 work = new long double[1+5*neqn]; 
 //file pointers has to be completed
 //file *sol, *sol_xy,

 //initial values 
 t=0.0;
 tend=1.48; 
 ans=ns+1.0;

 for(int i=0; i<ns; i++)
  for(int j=0; j<ns; j++)
   y[i*ns+j]=1.0;

 rtol=pow(10,-3);
 Atol=rtol;  //before it was small a now it is captal a because atol is a keywork in c++
 h=rtol;
 spcrad[0]=1;
 spcrad[1]=0;

 *idid=0;
 //need to start timing, 
  auto start = std::chrono::system_clock::now();
  ESERK(neqn, t, tend, h,   y, Atol,  spcrad, iwork,   work, idid);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  cout << "Time needed  " << elapsed_seconds.count() << endl;

 //ESERK(neqn, t, tend, dt, *g, atol, *spcrad,  *iwork, *work,  *idid_p)
 //need to end time


 //need to write the y values  to the file sol_serk5.out 
 //need to write the y values to the file sol_serk5_xy.out but only values at x=y
 //r is the distance
 
  //delete  idid;
  //delete [] iwork;
  //delete [] spcrad;
  delete [] y;
  //delete [] truey;
 // delete [] work;
return 0;
}

