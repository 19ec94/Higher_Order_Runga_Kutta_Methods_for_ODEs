#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>

inline double rho(int neqn,double t, double *y);
double ESERK(int neqn, double t, double tend, double dt, double *g, double tol,int *spcrad, int *iwork, double *work,int *idid_p);
double SERKrho(int neqn,double t, double *yn, double *fn,int *iwork,double hmax, double *work, double *sprad,int *idid);
double f(int neqn,double  x, double *y, double * f);

#endif /*FUNCTIONS_H*/
