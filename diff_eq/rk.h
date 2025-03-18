#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void rkf45(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);
void rk4(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);
void rk2(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);


void copy_vector(double * src, double * dest, int n);