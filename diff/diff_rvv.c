#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../gsl_rvv_math.h"
#include "gsl_rvv_diff.h"
#include <riscv-vector.h>

int
gsl_diff_backward (const gsl_function * f,
                   double x, double *result, double *abserr)
{

  int i, k;
  double h = GSL_SQRT_DBL_EPSILON;
  double a[3], d[3], a2;

  // i = 0, 1, 2
  // a = x - 2h, x - h, x
  // d = f(x - 2h), f(x - h), f(x)
  double koef[3] = {-2,-1,0};
  size_t vl = vsetvl_e64m2(3);
  vfloat64m2_t vec_a = vfmv_v_f_f64m2(x, vl);
  vfloat64m2_t vec_k = vle_v_f64m2(koef, vl);

  vec_a = vfmacc_vf_f64m2 (vec_a, h, vec_k, vl);
  vse_v_f64m2(a, vec_a, vl);

  

  for (i = 0; i < 3; i++)
    {
      d[i] = GSL_FN_EVAL (f, a[i]);
    }

  // k = 1, 2, 3

  for (k = 1; k < 4; k++)
    {
      // i = 0, 1
      // d[0] = (d[1]-d[0]) / (h)
      // d[1] = (d[2]-d[1]) / (h)
      
      // i = 0
      // d[0] = (d[1]-d[0]) / (2h) = (d[2] - 2d[1] + d[0])/2h^2
      for (i = 0; i < 3 - k; i++)
        {
          d[i] = (d[i + 1] - d[i]) / (a[i + k] - a[i]);
        }
    }

  
  //vlmax = vlen * lmul / sew
  // vl - длина регистра в элементах
  
  // sew = 64 (double) lmul = 2  vlen = 128 (const)
  vl = vsetvl_e64m2(3); 
  // становится равной 3, vlmax при этом 4.
  
  vfloat64m2_t vec_d = vle_v_f64m2(d, vl);
  vfloat64m1_t vec_sum = vfmv_v_f_f64m1(0.0, 1);
  vec_sum = vfredsum_vs_f64m2_f64m1(vec_sum, vec_d, vec_sum, vl);

  a2 = fabs(vfmv_f_s_f64m1_f64(vec_sum,1));

  if (a2 < 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      a2 = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  h = sqrt (GSL_SQRT_DBL_EPSILON / (2.0 * a2));

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      h = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  *result = (GSL_FN_EVAL (f, x) - GSL_FN_EVAL (f, x - h)) / h;
  *abserr = fabs (10.0 * a2 * h);

  return 0;
}

int
gsl_diff_forward (const gsl_function * f,
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
  double koef[3] = {0,1,2};
  size_t vl = vsetvl_e64m2(3);
  vfloat64m2_t vec_a = vfmv_v_f_f64m2(x, vl);
  vfloat64m2_t vec_k = vle_v_f64m2(koef, vl);

  vec_a = vfmacc_vf_f64m2 (vec_a, h, vec_k, vl);
  vse_v_f64m2(a, vec_a, vl);
  for (i = 0; i < 3; i++)
    {
      d[i] = GSL_FN_EVAL (f, a[i]);
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

  
  vl = vsetvl_e64m2(3); 
  
  vfloat64m2_t vec_d = vle_v_f64m2(d, vl);
  vfloat64m1_t vec_sum = vfmv_v_f_f64m1(0.0, 1);
  vec_sum = vfredsum_vs_f64m2_f64m1(vec_sum, vec_d, vec_sum, vl);

  a2 = fabs(vfmv_f_s_f64m1_f64(vec_sum,1));

  if (a2 < 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      a2 = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  h = sqrt (GSL_SQRT_DBL_EPSILON / (2.0 * a2));

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      h = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  *result = (GSL_FN_EVAL (f, x + h) - GSL_FN_EVAL (f, x)) / h;
  *abserr = fabs (10.0 * a2 * h);

  return 0;
}

int
gsl_diff_central (const gsl_function * f,
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
double koef[4] = {-2,-1,0,1};
  size_t vl = vsetvl_e64m2(4);
  vfloat64m2_t vec_a = vfmv_v_f_f64m2(x, vl);
  vfloat64m2_t vec_k = vle_v_f64m2(koef, vl);

  vec_a = vfmacc_vf_f64m2 (vec_a, h, vec_k, vl);
  vse_v_f64m2(a, vec_a, vl);
  for (i = 0; i < 4; i++)
    {
      a[i] = x + (i - 2.0) * h;
      d[i] = GSL_FN_EVAL (f, a[i]);
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

  size_t vlmax = vsetvlmax_e64m2();
  vl = vsetvl_e64m2(4);
  vfloat64m2_t vec_d = vle_v_f64m2(d, vl);
  vfloat64m1_t vec_sum = vfmv_v_f_f64m1(0.0, vl);
  vec_sum = vfredsum_vs_f64m2_f64m1(vec_sum, vec_d, vec_sum, vl);

  a3 = fabs(vfmv_f_s_f64m1_f64(vec_sum,1));

  
  //a3 = fabs (d[0] + d[1] + d[2] + d[3]);

  if (a3 < 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      a3 = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  h = pow (GSL_SQRT_DBL_EPSILON / (2.0 * a3), 1.0 / 3.0);

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    {
      h = 100.0 * GSL_SQRT_DBL_EPSILON;
    }

  *result = (GSL_FN_EVAL (f, x + h) - GSL_FN_EVAL (f, x - h)) / (2.0 * h);
  *abserr = fabs (100.0 * a3 * h * h);

  return 0;
}
