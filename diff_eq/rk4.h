#include <stdio.h>
#include <stdlib.h>

void runge_kutta(void (*f)(double, double*, double*), double t0, double * y0, double t_end, double h, int n);
