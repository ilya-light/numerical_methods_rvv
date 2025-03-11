#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "diff.h"

double f1 (double x)
{
  return exp (x);
}

double f2 (double x)
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

double f3 (double x)
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

int main ()
{
    double result_1, result_2, result_3 = 0;
    double abserr_1, abserr_2, abserr_3 = 0;

    for (size_t i = 0; i < 1000000; i++)
    {
      gsl_diff_backward(&f1, 1.5, &result_1, &abserr_1);
      gsl_diff_backward(&f2, 1.5, &result_2, &abserr_2);
      gsl_diff_backward(&f3, 1.5, &result_3, &abserr_3);
      gsl_diff_forward(&f1, 1.5, &result_1, &abserr_1);
      gsl_diff_forward(&f2, 1.5, &result_2, &abserr_2);
      gsl_diff_forward(&f3, 1.5, &result_3, &abserr_3);
      gsl_diff_central(&f1, 1.5, &result_1, &abserr_1);
      gsl_diff_central(&f2, 1.5, &result_2, &abserr_2);
      gsl_diff_central(&f3, 1.5, &result_3, &abserr_3);  
    }
}
