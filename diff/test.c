#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#undef GSL_DISABLE_DEPRECATED
#include "gsl_rvv_diff.h"

double
f1 (double x, void *params)
{
  return exp (x);
}

double
f2 (double x, void *params)
{
  if (x >= 0.0)
    {
      return 2 * x * x * x;
    }
  else
    {
      return 0.0;
    }
}

double
f3 (double x, void *params)
{
  if (x != 0.0)
    {
      return sin (1 / x);
    }
  else
    {
      return 0.0;
    }
}

int
main ()
{
  gsl_function F1, DF1, F2, DF2, F3, DF3, F4, DF4, F5, DF5, F6, DF6;

  F1.function = &f1;
  F2.function = &f2;
  F3.function = &f3;

  double result_1, result_2, result_3 = 0;
  double abserr_1, abserr_2, abserr_3 = 0;

  for (size_t i = 0; i < 1000000; i++)
  {
  gsl_diff_backward(&F1, 1.5, &result_1, &abserr_1);
  gsl_diff_backward(&F2, 1.5, &result_2, &abserr_2);
  gsl_diff_backward(&F3, 1.5, &result_3, &abserr_3);
  gsl_diff_forward(&F1, 1.5, &result_1, &abserr_1);
  gsl_diff_forward(&F2, 1.5, &result_2, &abserr_2);
  gsl_diff_forward(&F3, 1.5, &result_3, &abserr_3);
  gsl_diff_central(&F1, 1.5, &result_1, &abserr_1);
  gsl_diff_central(&F2, 1.5, &result_2, &abserr_2);
  gsl_diff_central(&F3, 1.5, &result_3, &abserr_3);  
  }
}
