
//#include "Constants1.h"
#include <iostream>
#include <cmath>
using namespace std;


//fcombustion(neqn,x,y,f)
//int main()
//{
  //neqn should be specified
long double fc(int neqn,long double x,long double *y,long double* f)
{
 int ns ,nsnsm1,i,ix,iy;
 long double ans,xx,yy; //y[neqn], f[neqn], x,;
 long double uxx,uyy,uij;
 long double d, alpha, delta, R,oneoverh2,par;
 
 /*Global veriables and their definitions */ 
 ns = sqrt(neqn);
 nsnsm1 =ns*(ns-1);
 oneoverh2 = (ns+1.0)*(ns+1.0);

 /*constants for inhomogenity */
  d=2.5;
  alpha = 1.0;
  delta =20.0;
  R=5.0;
  par=1.0;

//for(int i=0; i<neqn; i++){y[i]=1.0;} //Initialize the y
  for(int i=0; i<neqn; i++)
  {  
    if( (i%ns)  ==0)       {uxx = (-2.0*y[i]+2.0*y[i+1] )/ 3.0;}
    else if( (i+1)%ns ==0) {uxx = y[i-1]-2.0*y[i]+1.0;} 
    else                   {uxx = y[i-1]-2.0*y[i]+y[i+1];}
  
    if( i < ns)            {uyy = (-2.0*y[i]+2.0*y[i+ns] )/ 3.0;}
    else if(i >= nsnsm1 )  {uyy = y[i-ns]-2.0*y[i]+1.0;} 
    else                   {uyy = y[i-ns]-2.0*y[i]+y[i+ns];}


    uij = y[i];
    f[i]= d * oneoverh2*(uxx+uyy)+ (R/(alpha*delta) )*(1.0+alpha-uij)*exp(delta*(1.0-par/uij));  
  }
  /*for(int i=0; i<neqn ; i++)
  {
    cout<<"f["<<i<<"]="<<f[i]<<endl;
  }
  cout<<"end"<< endl;*/
return 0;

}

