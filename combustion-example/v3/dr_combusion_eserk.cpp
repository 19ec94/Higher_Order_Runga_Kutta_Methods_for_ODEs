
#include <iostream>
#include <cmath>
#include "functions.h"
#include <chrono>
#include <ctime>
#include <fstream>
#include "omp.h"
using namespace std;

template <int order>
inline void RTOL(double *rtol_p);

int main()
{
  int *idid = NULL, *iwork = NULL, *spcrad = NULL, index;
  const int ns = 400, neqn = ns * ns;
  double *y = NULL, error = 0.0;
  double *truey = NULL, ans = 0, xx = 0, yy = 0;
  double t = 0, tend = 0, r = 0;
  double h = 0;
  double rtol = 0, Atol = 0;
  double *work = NULL;

  idid = new int;
  iwork = new int[10];
  spcrad = new int[2];
  y = new double[neqn];
  truey = new double[neqn];
  work = new double[1 + 5 * neqn];
  FILE *sol, *sol_xy;

  sol = fopen("sol_serk5.txt", "w"); //NEEDS TO BE CHANGED according to ORDER
  sol_xy = fopen("sol_serk5_xy.txt", "w");

  t = 0.0;
  tend = 1.48;
  ans = ns + 1.0;

  for (int i = 0; i < neqn; i++)
  {
    y[i] = 1.0;
  }

  RTOL<ORDER>(&rtol);
  Atol = rtol;
  h = rtol;
  spcrad[0] = 1;
  spcrad[1] = 0;

  *idid = 0;
  /********************************
  Initialize the Integration
**********************************/
   double s = omp_get_wtime();
  //auto start = std::chrono::system_clock::now();
  ESERK(neqn, t, tend, h, y, Atol, spcrad, iwork, work, idid);
   double e = omp_get_wtime();
   cout <<"Time needed "<< e-s<<endl;
  //auto end = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end - start;
  //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  //cout << "Time needed  " << elapsed_seconds.count() << endl;

  /***********************************
  * Print solution to files
  * *********************************/
  for (int j = 0; j < neqn; j++)
  {
    fprintf(sol, "%d \t %.16f\n", j, y[j]);
  }
  fclose(sol);
  for (int j = 0; j < ns; j++)
  {
    xx = (j) / ans;
    yy = (j) / ans;
    r = sqrt((xx * xx + yy * yy));
    index = j * ns + j;
    fprintf(sol_xy, "%.17f \t %.16f\n", r, y[index]);
  }
  fclose(sol_xy);

  /*------------------------------------------------------
  Done.  Report some statistics.
------------------------------------------------------*/
  /* ----- print statistics -----*/
  cout << "Solution is tabulated in file sol_serk_xy.out" << endl;
  cout << "The value of IDID is " << *idid << endl;
  //cout<<"Max estimation of the spectral radius= "<<iwork[11]<<endl;
  //cout<<"Min estimation of the spectral radius= "<<iwork[12]<<endl;
  cout << "Max number of stages used=" << iwork[10] << endl;
  cout << "Number of f eval. for the spectr. radius=" << iwork[9] << endl;
  cout << "Number of f evaluations= " << iwork[5] << " steps= " << iwork[7] + iwork[8]
       << " accpt=" << iwork[7] << " rejct= " << iwork[8] << endl;

  /*------------------------------------------------------
  Done.  Compute the error and report some statistics.
------------------------------------------------------*/
  /*   
   std::ifstream infile("sol_exact.out");
   int i=0;
   while (infile >> truey[i]){}
   infile.close();
   
   for(int j=0; j<neqn; j++){error= max(error,abs(y[j]-truey[j]));}
   cout<<"error "<<error <<endl;
*/
  //delete idid;
  //delete[] iwork;
  //  delete [] spcrad;
  //    delete [] y;
  //delete [] truey;
  //  delete [] work;
  return 0;
}

template <int order>
inline void RTOL(double *rtol_p)
{
  if (order == 4)
  {
    *rtol_p = pow(10.0, -7);
  }
  else
  {
    *rtol_p = pow(10.0, -3);
  }
  return;
}
