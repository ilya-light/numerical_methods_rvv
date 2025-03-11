#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "diff.h"

int gsl_diff_backward (double (*f)(double),
                   double x, double *result, double *abserr)
{

  int i, k;
  double h = GSL_SQRT_DBL_EPSILON;
  double a[3], d[3], a2;

  for (i = 0; i < 3; i++)
    {
      a[i] = x + (i - 2.0) * h;
      d[i] = f(a[i]);
    }

  for (k = 1; k < 4; k++)
    {
      for (i = 0; i < 3 - k; i++)
        {
          d[i] = (d[i + 1] - d[i]) / (a[i + k] - a[i]);
        }
    }

  a2 = fabs (d[0] + d[1] + d[2]);

  if (a2 < 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      a2 = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  h = sqrt (GSL_SQRT_DBL_EPSILON / (2.0 * a2));

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      h = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  *result = (f(x) - f(x - h)) / h;
  *abserr = fabs (10.0 * a2 * h);

  return 0;
}

int gsl_diff_forward (double (*f)(double),
                  double x, double *result, double *abserr)
{
  /* Construct a divided difference table with a fairly large step
     size to get a very rough estimate of f''.  Use this to estimate
     the step size which will minimize the error in calculating f'. */

  int i, k;
  double h = GSL_SQRT_DBL_EPSILON;
  double a[3], d[3], a2;

  /* Algorithm based on description on pg. 204 of Conte and de Boor
     (CdB) - coefficients of Newton form of polynomial of degree 2. */

  for (i = 0; i < 3; i++)
    {
      a[i] = x + i * h;
      d[i] = f(a[i]);
    }

  for (k = 1; k < 4; k++)
    {
      for (i = 0; i < 3 - k; i++)
        {
          d[i] = (d[i + 1] - d[i]) / (a[i + k] - a[i]);
        }
    }

  /* Adapt procedure described on pg. 282 of CdB to find best value of
     step size. */

  a2 = fabs (d[0] + d[1] + d[2]);

  if (a2 < 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      a2 = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  h = sqrt (GSL_SQRT_DBL_EPSILON / (2.0 * a2));

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      h = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  *result = (f(x + h) - f(x)) / h;
  *abserr = fabs (10.0 * a2 * h);

  return 0;
}

int gsl_diff_central (double (*f)(double),
                  double x, double *result, double *abserr)
{
  /* Construct a divided difference table with a fairly large step
     size to get a very rough estimate of f'''.  Use this to estimate
     the step size which will minimize the error in calculating f'. */

  int i, k;
  double h = GSL_SQRT_DBL_EPSILON;
  double a[4], d[4], a3;

  /* Algorithm based on description on pg. 204 of Conte and de Boor
     (CdB) - coefficients of Newton form of polynomial of degree 3. */

  for (i = 0; i < 4; i++)
    {
      a[i] = x + (i - 2.0) * h;
      d[i] = f(a[i]);
    }

  for (k = 1; k < 5; k++)
    {
      for (i = 0; i < 4 - k; i++)
        {
          d[i] = (d[i + 1] - d[i]) / (a[i + k] - a[i]);
        }
    }

  /* Adapt procedure described on pg. 282 of CdB to find best
     value of step size. */

  a3 = fabs (d[0] + d[1] + d[2] + d[3]);

  if (a3 < 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      a3 = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  h = pow (GSL_SQRT_DBL_EPSILON / (2.0 * a3), 1.0 / 3.0);

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      h = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  *result = (f(x + h) - f(x - h)) / (2.0 * h);
  *abserr = fabs (100.0 * a3 * h * h);

  return 0;
}
