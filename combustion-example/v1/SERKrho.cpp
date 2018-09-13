#include <iostream>
#include <cmath>
#include <algorithm>
#include "functions.h"
using namespace std;

double SERKrho(int neqn, double t, double *yn, double *fn, int *iwork, double hmax, double *work, double *sprad, int *idid)
{
  const int itmax = 100;
  int index = 0, ptr5;
  double uround = 2.22e-16;
  double sqrtu = 0.0;
  double ynrm = 0.0;
  double sigma = 0.0, sigmal = 0.0, dynrm = 0.0, dfnrm = 0.0, vnrm = 0.0, small = 0.0;
  double *v = NULL, *fv = NULL;
  int nsteps = 0, nfesig = 0;
  sqrtu = sqrt((double)uround);
  nfesig = iwork[9];
  small = 1.0 / hmax;
  ptr5 = 4 * neqn;
  v = new double[neqn];
  fv = new double[neqn];

  nsteps = iwork[7] + iwork[8];
  for (int i = 0; i < neqn; i++)
  {
    v[i] = 0.0;
    fv[i] = 0.0;
  }

  if (nsteps == 0)
    for (int i = 0; i < neqn; i++)
    {
      v[i] = fn[i];
    }
  else
    for (int i = 0; i < neqn; i++)
    {
      v[i] = work[ptr5 + i];
    }

  ynrm = 0.0;
  vnrm = 0.0;

  for (int i = 0; i < neqn; i++)
  {
    ynrm = ynrm + pow(yn[i], 2);
    vnrm = vnrm + pow(v[i], 2);
  }
  ynrm = sqrt(ynrm);
  vnrm = sqrt(vnrm);
  if ((ynrm != 0.0) && (vnrm != 0.0))
  {
    dynrm = ynrm * sqrtu;
    for (int i = 0; i < neqn; i++)
    {
      v[i] = yn[i] + v[i] * (dynrm / vnrm);
    }
  }
  else if (ynrm != 0.0)
  {
    dynrm = ynrm * sqrtu;
    for (int i = 0; i < neqn; i++)
    {
      v[i] = yn[i] + yn[i] * sqrtu;
    }
  }
  else if (vnrm != 0.0)
  {
    dynrm = uround;
    for (int i = 0; i < neqn; i++)
    {
      v[i] = v[i] * (dynrm / vnrm);
    }
  }
  else
  {
    dynrm = uround;
    for (int i = 0; i < neqn; i++)
    {
      v[i] = dynrm;
    }
  }
  //Now iterate with a nonlinear power method
  sigma = 0.0;
  for (int iter = 1; iter <= itmax; iter++)
  {
    f(neqn, t, v, fv);
    nfesig = nfesig + 1;
    dfnrm = 0.0;
    for (int i = 0; i < neqn; i++)
    {
      dfnrm = dfnrm + pow((double)(fv[i] - fn[i]), 2);
    }
    dfnrm = sqrt((double)dfnrm);
    sigmal = sigma;
    sigma = dfnrm / dynrm;
    *sprad = 1.2 * sigma;
    if ((iter >= 2) && (fabs(sigma - sigmal) <= (max(sigma, small) * (0.01))))
    {
      for (int i = 0; i < neqn; i++)
      {
        work[i] = yn[i];
        work[neqn + i] = fn[i];
        work[2 * neqn + i] = v[i];
        work[3 * neqn + i] = fv[i];
        work[ptr5 + i] = v[i] - yn[i];
      }
      iwork[9] = nfesig;
      return 0;
    }
    if (dfnrm != 0.0)
    {
      for (int i = 0; i < neqn; i++)
      {
        v[i] = yn[i] + (fv[i] - fn[i]) * (dynrm / dfnrm);
      }
    }
    else
    {
      index = (iter % neqn);
      v[index] = yn[index] - (v[index] - yn[index]);
    }
  }
  iwork[9] = nfesig;
  *idid = 6;

  delete[] v;
  delete[] fv;

  return 0;
}
