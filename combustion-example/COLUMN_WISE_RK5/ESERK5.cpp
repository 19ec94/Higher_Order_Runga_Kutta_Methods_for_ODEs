#include <iostream>
#include <cmath>
#include <algorithm>
#include "Constants1.h"
#include <stdio.h>
using namespace std;
//old function parameters
//long double ESERK(neqn,t,tend,h,y,fcombusion, atol,spcrad, iwork,work,idid_p);
//New function parameters
//long double ESERK(neqn,t,tend,dt,g,f,tol,specrad,iwork,work,idid_p);

/*
int main()
{
//passed variables
//      integer          neqn, idid_p, iwork(10), spcrad(2)
//      long double precision t, tend, dt, tol, work(1+5*neqn), g(neqn)

int neqn=25;
 long double t =0;
 long double tend=1.48;
 long double dt; 
 long double tol;
 long double *g = new long double[neqn]; 
 long double *iwork = new long double[10];
 long double *work = new long double[1+5*neqn];
 long double *spcrad= new long double[2];
 int *idid_p;
 
*/



long double ESERK(int neqn,long double t,long double tend,long double dt,long double *g,long double tol,int *spcrad, int *iwork,long double *work,int *idid_p)
{
 int m,sec;
 int stage,stage_intern, max_stages;
 int total,accepted,rejected,rej_counter, feval;
 int start; //starting index of the coefficient array b
 int idid;
 long double eigmax,*eigmax_p=NULL, rho;  //previously it was eigmax
 //long double *b;  //it is in the header folder
 long double *sum =NULL, err, fac;
 long double *g_calc=NULL,*g_help=NULL, *g_save=NULL;
 long double *var_1=NULL, *var_2=NULL,*sc=NULL;
 long double *g_oth=NULL, *g_rej=NULL, r,pas, x0n, *y0n=NULL;
 long double dtmax,ans,hmax,tfail,al1;
 long double **g_work;
 long double **y;

 idid_p = &idid;
 //b= new long double[15484];
 eigmax_p=&eigmax;
 sum =new long double[neqn];
 g_calc= new long double[neqn];
 g_help=new long double[neqn];
 g_save= new long double[neqn];
 g_oth= new long double[neqn];
 var_1=new long double[neqn];
 var_2=new long double[neqn];
 sc = new long double[neqn];
 g_rej = new long double[neqn];
 y0n=new long double[neqn];
 g_work = new long double*[neqn];
 for(int i=0; i<neqn; i++) {g_work[i]= new long double[2001];}
 y=new long double *[neqn];
 for(int i=0; i<neqn; i++) {y[i]= new long double[5];}
 err=0.0;
 max_stages=0; 
 m=0; //need explanation but it is used as termiantion condition
 total=0; 
 accepted=0;
 rejected=0; 
 rej_counter=11; //why is it to 11 but not to 0
 feval=0;
 tfail = t;   //tfail equal to the starting time;
// cout<<"all the variables are initialized"<<endl;

 //for(int i=0; i<neqn; i++){g[i]=1.0;} //initialization has to be removed later
 //FILE *f;
 //f=fopen("g_help.txt","w");
std::cout<<fixed;

 fc(neqn, t,g,g_calc); //FUNCTION CALL
// for(int i =0; i<neqn; i++){cout<<g_calc[i];}
 cout<<endl;
 feval +=1;
 total +=1;
 //cout<<"first function call finished"<<endl;
 for(int i=0; i<neqn; i++)
 {
  sum[i]=0;
  g_help[i]=g[i];
  g_save[i]=g_calc[i];
  g_rej[i] = g[i];
  var_1[i]=g_calc[i];
//  cout<<"f["<<i<<"] "<<g_calc[i]<<endl;
 }
//cout<<"copying values finished,"<<endl;
 iwork[7]=0;
 iwork[8]=0;
 iwork[9]=0;
 hmax =abs(tend-t);
  
 if(spcrad[0]==0)
 {
  eigmax = rho; //rho(neqn,t,g);    the dummy function
 }
 else
 {
  SERKrho(neqn,t,g_help,g_calc,iwork,hmax,work,eigmax_p,idid_p);
  //cout<<"Internal rho "<<eigmax<<"\n";
  printf("Internal rho %.16Lf\n",eigmax);
  cout<<"Number of Func evaluations "<<iwork[9]<<"\n";
 }   
 
 dtmax = 0.98*2000.0*2000.0/eigmax;
 ans =sqrt(eigmax *dt/0.98);  //dt is the step size = h

 //cout<<"eigenmax is "<<eigmax<<endl;
 //cout<<"idid is "<<idid<<endl;
 //Need to set the stages, stage_inter and start i
 if(ans > 1800){
            stage = 2000;
            start = 13483;
            stage_intern = 200;}
      else if(ans > 1600){
            stage = 1800;
            start = 11682;
            stage_intern = 200;}
      else if(ans > 1400) {
            stage = 1600;
            start = 10081;
            stage_intern = 200;}
      else if(ans > 1200){
            stage = 1400;
            start = 8680;
            stage_intern = 200;}
      else if(ans > 1000) {
            stage = 1200;
            start = 7479;
            stage_intern = 200;}
      else if(ans > 900) {
            stage = 1000;
            start = 6478;
            stage_intern = 100;}
      else if(ans > 800) {
            stage = 900;
            start = 5577;
            stage_intern = 100;}
      else if(ans > 700) {
            stage = 800;
            start = 4776;
            stage_intern = 100;  }
      else if(ans > 600) {
            stage = 700;
            start = 4075;
            stage_intern = 100;}
      else if(ans > 500) {
            stage = 600;
            start = 3474;
            stage_intern = 100;  }
      else if(ans > 450) {
            stage = 500;
            start = 2973;
            stage_intern = 50  ;}
      else if(ans > 400) {
            stage = 450;
            start = 2522;
            stage_intern = 50  ;}
      else if(ans > 350) {
            stage = 400;
            start = 2121;
            stage_intern = 50  ;}
      else if(ans > 300) {
            stage = 350;
            start = 1770;
            stage_intern = 50  ;}
      else if(ans > 250) {
            stage = 300;
            start = 1469;
            stage_intern = 50  ;}
      else if(ans > 200) {
            stage = 250;
            start = 1218;
            stage_intern = 50  ;}
      else if(ans > 150) {
            stage = 200;
            start = 1017;
            stage_intern = 50  ;}
      else if(ans > 100) {
            stage = 150;
            start = 866;
            stage_intern = 50  ;}
      else if(ans > 90) {
            stage = 100;
            start = 765;
            stage_intern = 10  ;}
      else if(ans > 80) {
            stage = 90;
            start = 674;
            stage_intern = 10  ;}
      else if(ans > 70) {
            stage = 80;
            start = 593;
            stage_intern = 10  ;}
      else if(ans > 60) {
            stage = 70;
            start = 522;
            stage_intern = 10  ;}
      else if(ans > 50) {
            stage = 60;
            start = 461;
            stage_intern = 10  ;}
      else if(ans > 45) {
            stage = 50;
            start = 410;
            stage_intern = 5;  }
      else if(ans > 40) {
            stage = 45;
            start = 364;
            stage_intern = 5;  }
      else if(ans > 35) {
            stage = 40;
            start = 323;
            stage_intern = 5;  }
      else if(ans > 30) {
            stage = 35;
            start = 287;
            stage_intern = 5;  }
      else if(ans > 25) {
            stage = 30;
            start = 256;
            stage_intern = 5;  }
      else if(ans > 20) {
            stage = 25;
            start = 230;
            stage_intern = 5;  }
      else if(ans > 19) {
            stage = 20;
            start = 209;
            stage_intern = 2;  }
      else if(ans > 18) {
            stage = 19;
            start = 189;
            stage_intern = 2;  }
      else if(ans > 17) {
            stage = 18;
            start = 170;
            stage_intern = 2;  }
      else if(ans > 16) {
            stage = 17;
            start = 152;
            stage_intern = 2;  }
      else if(ans > 15) {
            stage = 16;
            start = 135;
            stage_intern = 2;  }
      else if(ans > 14) {
            stage = 15;
            start = 119;
            stage_intern = 2;  }
      else if(ans > 13) {
            stage = 14;
            start = 104;
            stage_intern = 2;  }
      else if(ans > 12) {
            stage = 13;
            start = 90;
            stage_intern = 2;  }
      else if(ans > 11) {
            stage = 12;
            start = 77;
            stage_intern = 2;  }
      else if(ans > 10) {
            stage = 11;
            start = 65;
            stage_intern = 2;  }
      else if(ans > 9) {
            stage = 10;
            start = 54;
            stage_intern = 2;  }
      else if(ans > 8) {
            stage = 9;
            start = 44;
            stage_intern = 2;  }
      else if(ans > 7) {
            stage = 8;
            start = 35;
            stage_intern = 2;  }
      else if(ans > 6) {
            stage = 7;
            start = 27;
            stage_intern = 2;  }
      else if(ans > 5) {
            stage = 6;
            start = 20;
            stage_intern = 2;  }
      else if(ans > 4) {
            stage = 5;
            start = 14;
            stage_intern = 2;  }
      else if(ans > 3) {
            stage = 4;
            start = 9;
            stage_intern = 2;  }
      else if(ans > 2){
            stage = 3;
            start = 5;
            stage_intern = 2;  }
      else if(ans > 1) {
            stage = 2;
            start = 2;
            stage_intern = 2;  }
      else {
            stage = 1;
            start = 0;
            stage_intern = 2;  }
      
 cout<<"stage intern "<<stage_intern<<endl;
 al1 = 1.0/(stage*stage*49.0/100.0);
 cout<<"stage "<<stage<<endl;
 //cout<<"al1 "<<al1<<endl;
 printf("al1 %.16Lf\n",al1);
 //runga kutta 
 while( ((t+dt) <= tend) && (m <= 1) && (total <= pow(10,7) ))
 {/*START-111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111*/
     cout<<"*****************"<<endl;     
  for(int i1=1; i1<=5; i1++) //[PROBLAMATIC]
  {/*START-222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222*/
    pas=dt/(i1);
    x0n=t;
  // std::cout<<std::scientific;
    //cout<<"pas  "<<pas<<endl;
    //cout<< "x0n "<<x0n<<endl;
    
    for(int l=0; l<neqn; l++){y0n[l]=g[l];}   //1-copy
    if(t>0)
    {
      fc(neqn,t,g,g_calc); //f is the fcombusion function
      feval=feval+1;
      for(int l=0; l<neqn; l++){g_save[l]=g_calc[l];} //2-copy
      //cout<<"gsave == gcalc"<<endl;
      //for(int l=0; l<neqn; l++){cout<<g_calc[l];}

    }  
    for(int j=1; j<=i1; j++)  //[problamatic]
    {/*START-333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333*/
      r=x0n;  //3-copy
      for(int l=0; l<neqn; l++){g_work[l][0]=y0n[l]; } //long double check , row and colum order
      //fprintf(f,"%f ",g_work[l][0]);}
        //  fprintf(f,"\n");
      //cout<<"g_work"<<endl;
      //for(int l=0;l<neqn; l++){cout<<g_work[l][0];}
      //cout<<endl;
      for(int i=1; i<=(stage/stage_intern); i++) //[PROBLAMATIC]
      {/* START-4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444*/
          //  cout<<"stage/Stage_intern "<<stage/stage_intern<<endl;
        if(j==1 && i==1)  //previously it was 1 now both changed to 0
        {
          //cout <<"i"<<i<<endl;
          for(int l=0; l<neqn; l++){g_calc[l]=g_save[l];} //fourth copy
          //for(int l=0; l<neqn; l++){cout<<"g_calc"<<g_calc[l];}
          //cout<<endl; //fourth copy
        }
        else
        {
          for(int l=0; l<neqn; l++){g_help[l]= g_work[l][1+(i-1)*stage_intern-1]; }//calculate and copy 5th copy
          
          fc(neqn,r,g_help,g_calc); //function call
          feval=feval+1;
         // for(int l=0; l<neqn; l++){cout<<"g_calc"<<g_calc[l];}
         // cout<<endl;
        }
        for(int l=0; l<neqn; l++)
        g_work[l][2+(i-1)*stage_intern-1]= g_work[l][1+(i-1)*stage_intern-1]+pas*al1*g_calc[l];   //calculate 6th
        //cout<<g_work[l][1+i*stage_intern]<<" ";
        
        //cout<<endl;
        //cout<<1+(i-1)*stage_intern-1<<endl;
        r = x0n + (1+stage_intern*stage_intern*(i-1))*pas*al1;  //r update (realted to 3-copy)
        //cout<<"r"<<r<<endl;
        for(int k=2; k<=stage_intern; k++)
        { 
          for(int l=0; l<neqn; l++){g_help[l]=g_work[l][1+k+(i-1)*stage_intern-1-1];}
         // cout<<g_help[l]<<endl;}  
         // cout<<1+k+(i)*stage_intern-1<<endl;
          fc(neqn,r,g_help,g_calc);  //function call
          feval = feval+1;
          //cout<<"feval "<<feval<<endl;
          for(int l=0; l<neqn; l++)
          {
            g_work[l][1+k+(i-1)*stage_intern-1] =2.0*g_work[l][1+k+(i-1)*
                stage_intern-1-1]-g_work[l][1+k+(i-1)*stage_intern-2-1]+2.0*pas*al1*g_calc[l];
           //fprintf(f,"%f ",g_work[l][1+k+(i-1)*stage_intern-1]);
        
          } //fprintf(f,"\n");  
          r=x0n+al1*(k*k+stage_intern*stage_intern*(i-1))*pas;
          //cout<<"r "<<r<<endl;
         }
      }/* START-4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444*/
    
      //if the stage is odd
      if((stage%stage_intern) ==1 )
      {
        for(int l=0; l<neqn; l++){g_help[l]=g_work[l][stage-1];}
        //for(int l=0; l<neqn; l++){cout<<g_help[l];}
       // cout<<"__________________________________"<<endl;
        fc(neqn,r,g_help,g_calc);  //function call
        feval=feval+1;
        //for(int l=0; l<neqn; l++){cout<<g_calc[l];}
        //cout<<"__________________________________"<<endl;
        for(int l=0; l<neqn; l++)
        {g_work[l][stage]=g_work[l][stage-1]+pas*al1*g_calc[l];}
        //cout<<g_work[l][stage];
        //cout<<endl;
      
      } 
      for(int k=0; k<stage+1; k++)
      {
        for(int l=0; l<neqn; l++){sum[l]=sum[l]+b[(start)+k]*g_work[l][k];} 

        //for(int l=0; l<neqn; l++){cout<<sum[l]<<" ";}
       //cout<<"__________________________________"<<endl;
      }
      
      for(int l=0; l<neqn; l++){y[l][i1-1]=sum[l]; y0n[l]=sum[l];sum[l]=0.0; } //ACHTUNG, Y SECOND INDEX
        x0n=x0n+pas;
     //   cout<<"x0n"<<x0n<<endl;
    }/*END-3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333*/
   } /*END-22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222*/
      err =0.0;
double c5 = 625.0 / 24.0, c4 = 1024.0 / 24.0, c3 = 486.0 / 24.0;
        double c2 = 64.0 / 24.0, c1 = 1.0 / 24.0;
      for(int l=0; l<neqn; l++)
      {
        g_oth[l]= g[l];
        /*g[l]= 625.0/24.0*y[l][4]-1024.0/24.0*y[l][3]+       //new
          486.0/24.0*y[l][2]-64.0/24.0*y[l][1]+
          1.0/24.0*y[l][0];*/
	g[l] = c5 * y[l][4] - c4 * y[l][3] + //new
                   c3 * y[l][2] - c2 * y[l][1] +
                   c1 * y[l][0];
        sc[l]=tol + tol * max(g[l],g_oth[l]);
        err  =err + pow( (( y[l][0]/24.0-4.0*y[l][1]/3.0+
                    27.0/4.0*y[l][2]-32.0/3.0*y[l][3]+
                    125.0/24.0*y[l][4] )*2.0/sc[l]),2); 
      //cout<<g_oth[l]<<" ";
      //cout<<g[l]<<" ";
      //cout<<sc[l]<< " ";
      }
      //cout<<endl;
      
      err = sqrt(err/neqn);
      
      t= t+dt;
      fc(neqn,t,g,g_calc);
      feval=feval+1;
      for(int i=0; i<neqn; i++){var_2[i]=g_calc[i];}
     // for(int i=0; i<neqn; i++){cout<<var_2[i]<<" ";}cout<<endl;
      
      //reject if error is too large
      if((isnan(err) || ((1.0/err) <=1.0)))
      {
        rejected=rejected+1;
        rej_counter=0;
        tfail =t;
        t =t-dt;
        fac = 0.8*pow((1.0/err),(1.0/5.0));
        //cout<<"fac "<<fac<<endl;
        fac = min(10.0L,fac);
        fac = max(fac,pow(10.0L,-4));
        //cout<<"fac "<<fac<<endl;
        dt = dt*fac;
        dt = min(dtmax, dt);
        //cout<<"dt "<<dt<<endl;
        //reset g
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
          fac = 0.8*pow((1.0/err),(1.0/5.0));
          if(t< tfail)
            sec=1;
          else
            sec=0;
      
          if((rej_counter <2) && sec==1)
          {
            rej_counter = rej_counter+1;
            fac = min(1.0L,fac);
          }
          else if ((rej_counter < 5) && (sec=1) )
          {
            rej_counter=rej_counter+1;
            fac=min(2.5L,fac); 
          }
          else
          {
            fac=min(10.0L,fac);
          }
          fac=max(fac,0.1L);
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
   dt = min( (0.98*2000.0*2000.0/eigmax),dt);
   ans = sqrt( eigmax*dt/0.98);
   
   printf("err = %.16Lf\n",err);
   printf("eigenmax %.16Lf\n",eigmax);
   printf("dt  %.16Lf\n",dt);
   printf("ans %.16Lf\n",ans);
   //cout<<"err "<<err<<endl;
   //cout<<"eigen max   "<<eigmax<<endl;
   //cout<<"dt "<<dt<<endl;
   //cout<<"ans "<<ans<<endl;
   
  // cout<<"feval "<<feval<<endl;
   //if else need to be added
    if(ans > 1800){
            stage = 2000;
            start = 13483;
            stage_intern = 200;  }
         else if(ans > 1600){
            stage = 1800;
            start = 11682;
            stage_intern = 200;  }
         else if(ans > 1400){
            stage = 1600;
            start = 10081;
            stage_intern = 200;}
         else if(ans > 1200){
            stage = 1400;
            start = 8680;
            stage_intern = 200;}
         else if(ans > 1000){
            stage = 1200;
            start = 7479;
            stage_intern = 200;}
         else if(ans > 900){
            stage = 1000;
            start = 6478;
            stage_intern = 100;}
         else if(ans > 800){
            stage = 900;
            start = 5577;
            stage_intern = 100;}
         else if(ans > 700){
            stage = 800;
            start = 4776;
            stage_intern = 100;}
         else if(ans > 600){
            stage = 700;
            start = 4075;
            stage_intern = 100;}
         else if(ans > 500){
            stage = 600;
            start = 3474;
            stage_intern = 100;}
         else if(ans > 450){
            stage = 500;
            start = 2973;
            stage_intern = 50;  }
         else if(ans > 400){
            stage = 450;
            start = 2522;
            stage_intern = 50 ; }
         else if(ans > 350){
            stage = 400;
            start = 2121;
            stage_intern = 50 ; }
         else if(ans > 300){
            stage = 350;
            start = 1770;
            stage_intern = 50;  }
         else if(ans > 250){
            stage = 300;
            start = 1469;
            stage_intern = 50;  }
         else if(ans > 200){
            stage = 250;
            start = 1218;
            stage_intern = 50 ; }
         else if(ans > 150){
            stage = 200;
            start = 1017;
            stage_intern = 50;  }
         else if(ans > 100){
            stage = 150;
            start = 866;
            stage_intern = 50 ; }
         else if(ans > 90){
            stage = 100;
            start = 765;
            stage_intern = 10;  }
         else if(ans > 80){
            stage = 90;
            start = 674;
            stage_intern = 10 ; }
         else if(ans > 70){
            stage = 80;
            start = 593;
            stage_intern = 10 ; }
         else if(ans > 60){
            stage = 70;
            start = 522;
            stage_intern = 10;  }
         else if(ans > 50){
            stage = 60;
            start = 461;
            stage_intern = 10 ; }
         else if(ans > 45){
            stage = 50;
            start = 410;
            stage_intern = 5;  }
         else if(ans > 40){
            stage = 45;
            start = 364;
            stage_intern = 5;  }
         else if(ans > 35){
            stage = 40;
            start = 323;
            stage_intern = 5;  }
         else if(ans > 30){
            stage = 35;
            start = 287;
            stage_intern = 5;  }
         else if(ans > 25){
            stage = 30;
            start = 256;
            stage_intern = 5;  }
         else if(ans > 20){
            stage = 25;
            start = 230;
            stage_intern = 5;  }
         else if(ans > 19){
            stage = 20;
            start = 209;
            stage_intern = 2;  }
         else if(ans > 18){
            stage = 19;
            start = 189;
            stage_intern = 2;  }
         else if(ans > 17){
            stage = 18;
            start = 170;
            stage_intern = 2;  }
         else if(ans > 16){
            stage = 17;
            start = 152;
            stage_intern = 2;  }
         else if(ans > 15){
            stage = 16;
            start = 135;
            stage_intern = 2;  }
         else if(ans > 14){
            stage = 15;
            start = 119;
            stage_intern = 2;  }
         else if(ans > 13){
            stage = 14;
            start = 104;
            stage_intern = 2;  }
         else if(ans > 12){
            stage = 13;
            start = 90;
            stage_intern = 2;  }
         else if(ans > 11){
            stage = 12;
            start = 77;
            stage_intern = 2;  }
         else if(ans > 10){
            stage = 11;
            start = 65;
            stage_intern = 2;  }
         else if(ans > 9){
            stage = 10;
            start = 54;
            stage_intern = 2;  }
         else if(ans > 8){
            stage = 9;
            start = 44;
            stage_intern = 2;  }
         else if(ans > 7){
            stage = 8;
            start = 35;
            stage_intern = 2;  }
         else if(ans > 6){
            stage = 7;
            start = 27;
            stage_intern = 2;  }
         else if(ans > 5){
            stage = 6;
            start = 20;
            stage_intern = 2;  }
         else if(ans > 4){
            stage = 5;
            start = 14;
            stage_intern = 2;  }
         else if(ans > 3){
            stage = 4;
            start = 9;
            stage_intern = 2;  }
         else if(ans > 2){
            stage = 3;
            start = 5;
            stage_intern = 2;  }
         else if(ans > 1){
            stage = 2;
            start = 2;
            stage_intern = 2;  }
         else{
            stage = 1;
            start = 0;
            stage_intern = 2;  }

      al1=1.0/(stage*stage*49.0/100.0);
      cout<<"stage "<<stage<<endl;
      cout<<"stage_intern "<<stage_intern<<endl;
      cout<<"start "<<start<<endl;
      //cout<<"al1 "<<al1<<endl;
      printf("al1 %.16Lf\n",al1);
      cout<<"*******************"<<endl;
 }/*END-1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111*/
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
