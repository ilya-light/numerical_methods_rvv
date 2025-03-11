#include "simple_plain.h"
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
    
    for (size_t i = 1; i < n; i++)
    {
        x = low + i * h;
        sum += f(x);
    }
    return sum * h;
}