#include "trapetional_integrate.h"
#include <math.h>
#include <riscv-vector.h>

double uniform_random() {
    return (double)rand() / (double)RAND_MAX;
}

double trapetional_integrate(const double (*f)(double),
                           const double low, const double up,
                           size_t n)
{
    double x, vol = 0;
    double h = (up - low) / n;
    double sum = 0.5 * (f(low) + f(up));
    
    n--;
    double * values = (double *)malloc(sizeof(double) * n);
    for (size_t i = 1; i <= n; i++)
    {
        values[i-1] = f(low + i * h);
    }
    size_t vlmax = vsetvlmax_e64m4();
    vfloat64m1_t vec_s = vfmv_v_f_f64m1(0, vlmax);
    vfloat64m4_t vec_values = vfmv_v_f_f64m4(0, vlmax);
    double * values_p = values;
    for (size_t vl; n > 0; n-=vl, values_p+=vl)
    {
        vl = vsetvl_e64m4(n);
        vec_values = vle_v_f64m4(values_p, vl);
        vec_s = vfredsum_vs_f64m4_f64m1(vec_s, vec_values, vec_s, vl);
    }
    
    free(values);
    return vfmv_f_s_f64m1_f64(vec_s,1) * h;

}