#include <stdlib.h>

#define POSINF (1.0 / 0.0)

double uniform_random();

double trapetional_integrate(const double (*f)(double),
                           const double low, const double up,
                           size_t n);
