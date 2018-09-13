#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


double f(int neqn,double x,double *y, double* f)
{
 int ns ,nsnsm1;
 double uxx,uyy,uij;
 double d, alpha, delta, R,oneoverh2,par;
 
 /*Global veriables and their definitions */ 
 ns = sqrt((double)neqn);
 nsnsm1 =ns*(ns-1);
 oneoverh2 = (ns+1.0L)*(ns+1.0L);
  d=2.5L;
  alpha = 1.0L;
  delta =20.0L;
  R=5.0L;
  par=1.0L;

  for(int i=0; i<neqn; i++)
  {  
    if( (i%ns)  ==0)       {uxx = (-2.0*y[i]+2.0*y[i+1] )/ 3.0;}
    else if( (i+1)%ns ==0) {uxx = y[i-1]-2.0*y[i]+1.0;} 
    else                   {uxx = y[i-1]-2.0*y[i]+y[i+1];}
  
    if( i < ns)            {uyy = (-2.0*y[i]+2.0*y[i+ns] )/ 3.0;}
    else if(i >= nsnsm1 )  {uyy = y[i-ns]-2.0*y[i]+1.0;} 
    else                   {uyy = y[i-ns]-2.0*y[i]+y[i+ns];}


    f[i]= d * oneoverh2*(uxx+uyy)+ (R/(alpha*delta) )*(1.0+alpha-y[i])*exp(delta*(1.0-par/y[i]));  
  }
  
return 0;

}
