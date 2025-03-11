#define GSL_SQRT_DBL_EPSILON 1.490116119384777e-8

int gsl_diff_central (double (*f)(double),
                      double x,
                      double *result, double *abserr);

int gsl_diff_backward (double (*f)(double),
                       double x,
                       double *result, double *abserr);

int gsl_diff_forward (double (*f)(double),
                      double x,
                      double *result, double *abserr);

