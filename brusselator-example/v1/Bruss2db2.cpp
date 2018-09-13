#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


double f(int neqn,double x,double *y, double* f)
//int main()
{
 int ns ,nsnsm1, nssq;
 double ans,xx,yy;
 double uleft,vleft,uright,vright,ulow,vlow;
 double uup,vup,uij,vij,brussa,brussb,alf;

/*
 int neqn = 2*2*2;
 double *y=new double [neqn];
 double *f=new double [neqn];
 for(int i=0; i<neqn; i++){y[i]=1.0;}
*/
 /*Global veriables and their definitions */ 
 
 nssq = neqn/2;
 ns = sqrt((double)nssq);
 nsnsm1 =ns*(ns-1);

 /*constants for inhomogenity*/
 brussa =1.3;
 brussb=2.0e6;
 alf=1.0e-1;
 
  // ----- big loop -----
      for(int i=1; i<=nssq; i++)
      {
// ----- left neighbour -----
         if((i% ns) == 1){
            uleft=y [(i+ns-1)*2-2];
            vleft=y [(i+ns-1)*2-1];
         }
         else
         {
            uleft=y [(i-1)*2-2 ];
            vleft=y [(i-1)*2-1];
         }
// ----- right neighbour -----
         if((i%ns) == 0)
         {
            uright=y[(i-ns+1)*2-2];
            vright=y[(i-ns+1)*2-1];
         }
         else
         {
            uright=y[(i+1)*2-2];
            vright=y[(i+1)*2-1];
         }
// ----- lower neighbour -----
         if(i <= ns)
         {
            ulow=y[(i+nsnsm1)*2-2];
            vlow=y[(i+nsnsm1)*2-1];
         }
         else
         {
            ulow=y[(i-ns)*2-2];
            vlow=y[(i-ns)*2-1];
         }
// ----- upper neighbour -----
         if(i > nsnsm1)
         {
            uup=y[(i-nsnsm1)*2-2];
            vup=y[(i-nsnsm1)*2-1];
         }
         else
         {
            uup=y[(i+ns)*2-2];
            vup=y[(i+ns)*2-1];
         }
// ----- the derivative -----
         uij=y [i*2-2 ];
         vij=y [i*2-1];
         f[i*2-2]=alf*nssq*(uleft+uright+ulow+uup-4.0*uij)+brussa+uij*uij*vij-(brussb+1.0)*uij ;
         f[i*2-1]=alf*nssq*(vleft+vright+vlow+vup-4.0*vij)+brussb*uij - uij*uij*vij ;
         
      }
     // std::cout<<fixed;
     // for(int i=0; i<neqn; i++) cout<<f[i]<<endl;

return 0;

}
