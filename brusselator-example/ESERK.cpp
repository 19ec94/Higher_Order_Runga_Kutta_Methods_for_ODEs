#include <iostream>
#include <cmath>
#include <algorithm>
#include "functions.h"
#if (ORDER==4)
#include "Constants4.h"
#else
#include "Constants5.h" 
#endif
#include <stdio.h>
using namespace std;

template<int order>
inline void ORDER_SELECTION(int *local_order);
template<int order>
inline void DTMAX(double *dtmax_p, double *eigmax_p);
template<int order>
inline void ANS(double *ans_p,double *dt_p, double *eigmax_p);
template<int order>
inline void STAGE_SELECTION(double *ans_p, int *stage_p, int *start_p, int *stage_intern_p);
template<int order>
inline void AL1(double *al1_p, int *stage_p );
template<int order>
inline void DT(double *eigmax_p, double *dt_p);
template<int order>
inline void ERROR_CALCULATION(int neqn, double *g_oth, double *g, double **y ,double *sc,double tol, double *err_p  );
template<int order>
inline void FAC(double *fac_p,double *err_p);
template<int order>
inline void ODD_STAGE(int stage,int stage_intern,int neqn,double *g_help, double **g_work,double *g_calc, double r, int *feval_p, double pas, double al1);


double ESERK(int neqn,double t,double tend,double dt, double *g, double tol,int *spcrad, int *iwork, double *work,int *idid_p)
{
 int m,sec;
 int stage,stage_intern, max_stages;
 int total,accepted,rejected,rej_counter, feval;
 int start; //starting index of the coefficient array b
 int idid;
 double eigmax,*eigmax_p=NULL, rho=0;  //previously it was eigmax
 double *sum =NULL;
 double err, fac;
 double *g_calc=NULL,*g_help=NULL, *g_save=NULL;
 double *var_1=NULL, *var_2=NULL,*sc=NULL;
 double *g_oth=NULL, *g_rej=NULL;
 double  r,pas, x0n;
 double *y0n=NULL;
 double dtmax;
 double ans;
 double hmax;
 double tfail,al1;
 double **g_work;
 double **y;
 int highest_order=0; 

 idid_p = &idid;
 eigmax_p=&eigmax;
 sum =new double[neqn];
 g_calc= new double[neqn];
 g_help=new double[neqn];
 g_save= new double[neqn];
 g_oth= new double[neqn];
 var_1=new double[neqn];
 var_2=new double[neqn];
 sc = new double[neqn];
 g_rej = new double[neqn];
 y0n=new double[neqn];
 g_work = new double*[2001];
 for(int i=0; i<2001; i++) {g_work[i]= new double[neqn];}
 
 ORDER_SELECTION<ORDER>(&highest_order);

 y=new double *[highest_order];
 for(int i=0; i<highest_order; i++) {y[i]= new double[neqn];}
 
 
 err=0.0;
 max_stages=0; 
 m=0; //need explanation but it is used as termiantion condition
 total=0; 
 accepted=0;
 rejected=0; 
 rej_counter=11; 
 feval=0;
 tfail = t;   //tfail equal to the starting time;
 
 f(neqn, t,g,g_calc); //FUNCTION CALL
 feval =feval+1;
 total =total+1;
 for(int i=0; i<neqn; i++)
 {
  sum[i]=0.0;
  g_help[i]=g[i];
  g_save[i]=g_calc[i];
  g_rej[i] = g[i];
  var_1[i]=g_calc[i];
 }
 iwork[7]=0;
 iwork[8]=0;
 iwork[9]=0;
 
 hmax =fabs(tend-t);
  
   

 if(spcrad[0]==0)
 {
  eigmax = rho; //rho(neqn,t,g);    the dummy function
 }
 else
 {
  SERKrho(neqn,t,g_help,g_calc,iwork,hmax,work,eigmax_p,idid_p);
  //cout<<"Internal rho "<<eigmax<<"\n";
  printf("Internal rho %.16f\n",eigmax);
  cout<<"Number of Func evaluations "<<iwork[9]<<"\n";
 }   
 

 DTMAX<ORDER>(&dtmax, &eigmax);
 ANS<ORDER>(&ans,&dt, &eigmax);
 STAGE_SELECTION<ORDER>(&ans, &stage, &start, &stage_intern);
 AL1<ORDER>(&al1, &stage );


 

 while( ((t+dt) <= tend) && (m <= 1) && (total <= pow(10,7) ))
 {
  for(int i1=1; i1<=highest_order; i1++) 
  {
    pas=dt/(double)i1;
    x0n=t;
    for(int l=0; l<neqn; l++){y0n[l]=g[l];}   //1-copy
    if(t>0)
    {
      f(neqn,t,g,g_calc); //f is the fcombusion function
      feval=feval+1;
      for(int l=0; l<neqn; l++){g_save[l]=g_calc[l];} //2-copy
      
    }  
    for(int j=1; j<=i1; j++) 
    {
       r=x0n;  //3-copy
      for(int l=0; l<neqn; l++){g_work[0][l]=y0n[l]; } 
      for(int i=1; i<=(stage/stage_intern); i++) 
      {
        if(j==1 && i==1) 
        {
          for(int l=0; l<neqn; l++){g_calc[l]=g_save[l];} //fourth copy
        }
        else
        {
          for(int l=0; l<neqn; l++){g_help[l]= g_work[1+(i-1)*stage_intern-1][l]; }//calculate and copy 5th copy
          f(neqn,r,g_help,g_calc); //function call
          feval=feval+1;
         }

        for(int l=0; l<neqn; l++){
        g_work[2+(i-1)*stage_intern-1][l]= g_work[1+(i-1)*stage_intern-1][l]+pas*al1*g_calc[l];   //calculate 6th
        }
        
       r = x0n + (1+stage_intern*stage_intern*(i-1))*pas*al1;  //r update (realted to 3-copy)
        for(int k=2; k<=stage_intern; k++)
        { 
          for(int l=0; l<neqn; l++){g_help[l]=g_work[1+k+(i-1)*stage_intern-1-1][l];}
          f(neqn,r,g_help,g_calc);  //function call
          feval = feval+1;
          for(int l=0; l<neqn; l++)
          {
            g_work[1+k+(i-1)*stage_intern-1][l] =2.0*g_work[1+k+(i-1)*
                stage_intern-1-1][l]-g_work[1+k+(i-1)*stage_intern-2-1][l]+2.0*pas*al1*g_calc[l];
          } 
          r=x0n+al1*(k*k+stage_intern*stage_intern*(i-1))*pas;
         }
      }
     
     ODD_STAGE<ORDER>(stage,stage_intern,neqn,g_help,g_work,g_calc, r, &feval, pas, al1);

     
      for(int k=0; k<stage+1; k++)
      {
        for(int l=0; l<neqn; l++){sum[l]=sum[l]+b[(start)+k]*g_work[k][l];} 
       }  


      for(int l=0; l<neqn; l++){y[i1-1][l]=sum[l]; y0n[l]=sum[l];sum[l]=0.0; } 
        x0n=x0n+pas;
      }
   } 
      err =0.0;
    ERROR_CALCULATION<ORDER>(neqn, g_oth, g, y ,sc,tol, &err  );
      err = sqrt(err/neqn);
      t= t+dt;
      f(neqn,t,g,g_calc);
      feval=feval+1;
      for(int i=0; i<neqn; i++){var_2[i]=g_calc[i];}
      
      //reject if error is too large
      if((isnan(err) || ((1.0/err) <=1.0)))
      {
        rejected=rejected+1;
        rej_counter=0;
        tfail =t;
        t =t-dt;
        FAC<ORDER>(&fac,&err);
        fac = min(10.0,fac);
        fac = max(fac,pow(10.0,-4));
        dt = dt*fac;
        dt = min(dtmax, dt);
        for(int i=0; i<neqn; i++)
        {
          g[i]=g_rej[i];
          g_help[i]=g[i];
          g_calc[i]=var_1[i];
        }
      }
      else
      {
        accepted=accepted+1;;
        if(total>1)
        {
          FAC<ORDER>(&fac, &err);
          if(t< tfail)
            sec=1;
          else
            sec=0;
      
          if((rej_counter <2) && sec==1)
          {
            rej_counter = rej_counter+1;
            fac = min(1.0,fac);
          }
          else if ((rej_counter < 5) && (sec=1) )
          {
            rej_counter=rej_counter+1;
            fac=min(2.5,fac); 
          }
          else
          {
            fac=min(10.0,fac);
          }
          fac=max(fac,0.1);
          dt=dt*fac;
          dt=min(dtmax,dt);
          }
    //accept step just calculated
        for(int i=0; i<neqn; i++)
        {
          g_rej[i]=g[i];
          g_help[i]=g[i];
          g_calc[i]=var_2[i];
          var_1[i]=g_calc[i];
        }
        if(stage > max_stages)
           max_stages=stage;
      }
      total=total+1;
      //shorten last step if it is to big
      if((t+dt) > tend)
      {
        dt=tend-t;
        m=m+1;
      }

  /* choose number of stages to use for next step [[[[[ACHTUNG - INDEX IS 0 AND 1 , BUT NOT 1 AND 2]]]]
     if h_new too bigd0, make it smaller
     spcrad(1) = 0  -> subroutine rho is providing an upper bound on the spectralradius
               = 1  -> compute an upper bound on the spectralradius using a nonlinear power method
     spcrad(2) = 0  ->  The Jacobian may not be constant.
               = 1  ->  The Jacobian is constant. */
   if(spcrad[1]==0)
   {
     if(spcrad[0]==0)
     {
       eigmax= rho; //rho(neqn,t,g);
     }
     else
     {
       SERKrho(neqn,t,g_help,g_calc,iwork,hmax,work, eigmax_p,idid_p);
     }
   }

  
    
 DT<ORDER>( &eigmax, &dt);
 ANS<ORDER>(&ans,&dt, &eigmax);
 STAGE_SELECTION<ORDER>(&ans, &stage, &start, &stage_intern);
 AL1<ORDER>(&al1, &stage );

      printf("eigenmax= %.23f \n",eigmax);
      //printf("dt = %.16e \n",dt);
      //cout<<"stage_intern "<<stage_intern<<endl;
      //cout<<"stage "<<stage<<endl;
      //cout<<"al1 "<<al1<<endl;
      //printf("al1 %.16e\n",al1);
      //cout<<"*******************"<<endl;
 }
 
  *idid_p =1;
  iwork[5]=feval;
  iwork[7]=accepted;
  iwork[8]=rejected;
  iwork[10]=max_stages;
  cout<<feval<<endl;
  delete [] sum;
  delete [] g_calc;
  delete [] g_help;
  delete [] g_save;
  delete [] var_1;
  delete [] var_2;
  delete [] sc;
  delete [] g_oth;
  delete [] g_rej;
  delete [] y0n;
 return 0;
}

template<int order>
inline void ORDER_SELECTION(int *local_order)
 {
       if(order ==4) {*local_order =4;}
       else          {*local_order =5;}
       return ;
 }

template<int order>
inline void DTMAX(double *dtmax_p, double *eigmax_p)
{
      if(order ==4) {*dtmax_p = 2000.0*2000.0/ (*eigmax_p);}
      else          {*dtmax_p = 0.98*2000.0*2000.0/(*eigmax_p);}
 
 return ;
}

template<int order>
inline void ANS(double *ans_p,double *dt_p, double *eigmax_p)
{
     if(order ==4) { *ans_p =sqrt( (*eigmax_p) * (*dt_p)); }
     else          {*ans_p =sqrt( (*eigmax_p) * (*dt_p)/0.98); }  //dt is the step size = h
    
    return ;
}

template<int order>
inline void STAGE_SELECTION(double *ans_p, int *stage_p, int *start_p, int *stage_intern_p)
{
if(order ==4)
{
         if((*ans_p) > 1800) {
            *stage_p = 2000;
            *start_p = 13265;
            *stage_intern_p = 200  ;}
         else if((*ans_p) > 1600) {
            *stage_p = 1800       ;
            *start_p = 11464      ;
            *stage_intern_p = 200  ;}
         else if((*ans_p) > 1400) {
            *stage_p = 1600       ;
            *start_p = 9863      ;
            *stage_intern_p = 200  ;}
         else if((*ans_p) > 1200) {
            *stage_p = 1400       ;
            *start_p = 8462      ;
            *stage_intern_p = 200  ;}
         else if((*ans_p) > 1000) {
            *stage_p = 1200       ;
            *start_p = 7261      ;
            *stage_intern_p = 200  ;}
         else if((*ans_p) > 900) {
            *stage_p = 1000       ;
            *start_p = 6260      ;
            *stage_intern_p = 100  ;}
         else if((*ans_p) > 800) {
            *stage_p = 900       ;
            *start_p = 5359      ;
            *stage_intern_p = 100  ;}
         else if((*ans_p) > 700) {
            *stage_p = 800       ;
            *start_p = 4558      ;
            *stage_intern_p = 100  ;}
         else if((*ans_p) > 600) {
            *stage_p = 700       ;
            *start_p = 3857      ;
            *stage_intern_p = 100  ;}
         else if((*ans_p) > 500) {
            *stage_p = 600       ;
            *start_p = 3256      ;
            *stage_intern_p = 100  ;}
         else if((*ans_p) > 450) {
            *stage_p = 500       ;
            *start_p = 2755      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 400) {
            *stage_p = 450       ;
            *start_p = 2304      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 350) {
            *stage_p = 400       ;
            *start_p = 1903      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 300) {
            *stage_p = 350       ;
            *start_p = 1552      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 250) {
            *stage_p = 300       ;
            *start_p = 1251      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 200) {
            *stage_p = 250       ;
            *start_p = 1000      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 150) {
            *stage_p = 200       ;
            *start_p = 799      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 100) {
            *stage_p = 150       ;
            *start_p = 648      ;
            *stage_intern_p = 25  ;}
         else if((*ans_p) > 90) {
            *stage_p = 100       ;
            *start_p = 547      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 80) {
            *stage_p = 90       ;
            *start_p = 456      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 70) {
            *stage_p = 80       ;
            *start_p = 375      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 60) {
            *stage_p = 70       ;
            *start_p = 304      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 50) {
            *stage_p = 60       ;
            *start_p = 243      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 40) {
            *stage_p = 50       ;
            *start_p = 192      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 30) {
            *stage_p = 40       ;
            *start_p = 151      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 20) {
            *stage_p = 30       ;
            *start_p = 120      ;
            *stage_intern_p = 10  ;}
         else if((*ans_p) > 18) {
            *stage_p = 20       ;
            *start_p = 99      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 16) {
            *stage_p = 18       ;
            *start_p = 80      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 14) {
            *stage_p = 16       ;
            *start_p = 63      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 12) {
            *stage_p = 14       ;
            *start_p = 48      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 10) {
            *stage_p = 12       ;
            *start_p = 35      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 8) {
            *stage_p = 10       ;
            *start_p = 24      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 6) {
            *stage_p = 8       ;
            *start_p = 15      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 4) {
            *stage_p = 6       ;
            *start_p = 8      ;
            *stage_intern_p = 2  ;}
         else if((*ans_p) > 2) {
            *stage_p = 4       ;
            *start_p = 3      ;
            *stage_intern_p = 2  ;}
         else{
            *stage_p = 2;
            *start_p = 0;
            *stage_intern_p = 2;}
}
else 
{
if( ((*ans_p)) > 1800){
            *stage_p = 2000;
            *start_p = 13483;
            *stage_intern_p = 200;}
      else if((*ans_p) > 1600){
            *stage_p = 1800;
            *start_p = 11682;
            *stage_intern_p = 200;}
      else if((*ans_p) > 1400) {
            *stage_p = 1600;
            *start_p = 10081;
            *stage_intern_p = 200;}
      else if( (*ans_p) > 1200){
            *stage_p = 1400;
            *start_p = 8680;
            *stage_intern_p = 200;}
      else if((*ans_p) > 1000) {
            *stage_p = 1200;
            *start_p = 7479;
            *stage_intern_p = 200;}
      else if((*ans_p) > 900) {
            *stage_p = 1000;
            *start_p = 6478;
            *stage_intern_p = 100;}
      else if((*ans_p) > 800) {
            *stage_p = 900;
            *start_p = 5577;
            *stage_intern_p = 100;}
      else if((*ans_p) > 700) {
            *stage_p = 800;
            *start_p = 4776;
            *stage_intern_p = 100;  }
      else if((*ans_p) > 600) {
            *stage_p = 700;
            *start_p = 4075;
            *stage_intern_p = 100;}
      else if((*ans_p) > 500) {
            *stage_p = 600;
            *start_p = 3474;
            *stage_intern_p = 100;  }
      else if((*ans_p) > 450) {
            *stage_p = 500;
            *start_p = 2973;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 400) {
            *stage_p = 450;
            *start_p = 2522;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 350) {
            *stage_p = 400;
            *start_p = 2121;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 300) {
            *stage_p = 350;
            *start_p = 1770;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 250) {
            *stage_p = 300;
            *start_p = 1469;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 200) {
            *stage_p = 250;
            *start_p = 1218;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 150) {
            *stage_p = 200;
            *start_p = 1017;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 100) {
            *stage_p = 150;
            *start_p = 866;
            *stage_intern_p = 50  ;}
      else if((*ans_p) > 90) {
            *stage_p = 100;
            *start_p = 765;
            *stage_intern_p = 10  ;}
      else if((*ans_p) > 80) {
            *stage_p = 90;
            *start_p = 674;
            *stage_intern_p = 10  ;}
      else if((*ans_p) > 70) {
            *stage_p = 80;
            *start_p = 593;
            *stage_intern_p = 10  ;}
      else if((*ans_p) > 60) {
            *stage_p = 70;
            *start_p = 522;
            *stage_intern_p = 10  ;}
      else if((*ans_p) > 50) {
            *stage_p = 60;
            *start_p = 461;
            *stage_intern_p = 10  ;}
      else if((*ans_p) > 45) {
            *stage_p = 50;
            *start_p = 410;
            *stage_intern_p = 5;  }
      else if((*ans_p) > 40) {
            *stage_p = 45;
            *start_p = 364;
            *stage_intern_p = 5;  }
      else if((*ans_p) > 35) {
            *stage_p = 40;
            *start_p = 323;
            *stage_intern_p = 5;  }
      else if((*ans_p) > 30) {
            *stage_p = 35;
            *start_p = 287;
            *stage_intern_p = 5;  }
      else if((*ans_p) > 25) {
            *stage_p = 30;
            *start_p = 256;
            *stage_intern_p = 5;  }
      else if((*ans_p) > 20) {
            *stage_p = 25;
            *start_p = 230;
            *stage_intern_p = 5;  }
      else if((*ans_p) > 19) {
            *stage_p = 20;
            *start_p = 209;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 18) {
            *stage_p = 19;
            *start_p = 189;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 17) {
            *stage_p = 18;
            *start_p = 170;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 16) {
            *stage_p = 17;
            *start_p = 152;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 15) {
            *stage_p = 16;
            *start_p = 135;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 14) {
            *stage_p = 15;
            *start_p = 119;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 13) {
            *stage_p = 14;
            *start_p = 104;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 12) {
            *stage_p = 13;
            *start_p = 90;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 11) {
            *stage_p = 12;
            *start_p = 77;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 10) {
            *stage_p = 11;
            *start_p = 65;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 9) {
            *stage_p = 10;
            *start_p = 54;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 8) {
            *stage_p = 9;
            *start_p = 44;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 7) {
            *stage_p = 8;
            *start_p = 35;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 6) {
            *stage_p = 7;
            *start_p = 27;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 5) {
            *stage_p = 6;
            *start_p = 20;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 4) {
            *stage_p = 5;
            *start_p = 14;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 3) {
            *stage_p = 4;
            *start_p = 9;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 2){
            *stage_p = 3;
            *start_p = 5;
            *stage_intern_p = 2;  }
      else if((*ans_p) > 1) {
            *stage_p = 2;
            *start_p = 2;
            *stage_intern_p = 2;  }
      else {
            *stage_p = 1;
            *start_p = 0;
            *stage_intern_p = 2;  }
}
 
 
      return;
}

template<int order>
 inline void AL1(double *al1_p, int *stage_p )
 {   
       if(order ==4) { *al1_p = 1.0/( (*stage_p) * (*stage_p)  /2.0);}
       else          { *al1_p = 1.0/( (*stage_p) * (*stage_p) *49.0/100.0);}
      
       return;
 }

template<int order>
   inline void DT(double *eigmax_p, double *dt_p)
   {
      if (order ==4){ *dt_p  = min( ( (2000.0*2000.0)/ (*eigmax_p)  ), (*dt_p) );}
      else          { *dt_p  = min( ( (0.98*2000.0*2000.0)/ (*eigmax_p) ),(*dt_p) ); }
         
         return ;
   }


template<int order>
   inline void ERROR_CALCULATION(int neqn, double *g_oth, double *g, double **y ,double *sc,double tol, double *err_p  )
{    

   if(order ==4)
   {
      for(int l=0; l<neqn; l++)
      {
      g_oth[l]= g[l];
      g[l]= 64.0/6.0*y[3][l]-81.0/6.0*y[2][l]
            +24.0/6.0*y[1][l]-1.0/6.0*y[0][l];
      sc[l]=tol + tol * max(g[l],g_oth[l]);
      *err_p  = (*err_p) + pow( (( -y[0][l]/6.0+2.0*y[1][l]-9.0/2.0*y[2][l]
                  +8.0/3.0*y[3][l])*2.0/sc[l]),2); 
      }
    }
  
   else 
    {
      double c5=625.0/24.0,c4=1024.0/24.0,c3=486.0/24.0;
      double c2=64.0/24.0,c1=1.0/24.0;
      for(int l=0; l<neqn; l++)
      {
        g_oth[l]= g[l];
        g[l]= c5*y[4][l]-c4*y[3][l]+       //new
          c3*y[2][l]-c2*y[1][l]+
          c1*y[0][l];
        sc[l]=tol + tol * max(g[l],g_oth[l]);
        *err_p  = (*err_p) + pow( (( y[0][l]/24.0-4.0*y[1][l]/3.0+
                    27.0/4.0*y[2][l]-32.0/3.0*y[3][l]+
                    125.0/24.0*y[4][l] )*2.0/sc[l]),2); 
      }
    }
    return ;

}
template<int order>
inline void FAC(double *fac_p,double *err_p)
{     
      if(order ==4){*fac_p = 0.8*pow((1.0/ (*err_p) ),(0.25));}
      else         {*fac_p = (8.0/10.0)*pow((1.0/ (*err_p) ),(1.0/5.0));}
      
      return ;
}

template<int order>
 inline void ODD_STAGE(int stage,int stage_intern,int neqn,double *g_help, double **g_work,double *g_calc, double r, int *feval_p, double pas, double al1)
{    if(order ==5)
     {
      if((stage%stage_intern) ==1 )
      {
        for(int l=0; l<neqn; l++){g_help[l]=g_work[stage-1][l];}
        f(neqn,r,g_help,g_calc);  //function call
        *feval_p = (*feval_p )+1;
        for(int l=0; l<neqn; l++)
        {g_work[stage][l]=g_work[stage-1][l]+pas*al1*g_calc[l];}
      }
     }   
}
