#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void rkf45(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);
void rk4(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);
void rk2(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);

void rkf45_rvv(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);
void rk4_rvv(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);
void rk2_rvv(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n);


void copy_vector(double * src, double * dest, int n);
void copy_vector_rvv(double * src, double * dest, int n);

void rkf45_rvv_1(int, double, double *, double *, double *);
void rkf45_rvv_2(int, double, double *, double *, double *, double *);
void rkf45_rvv_3(int, double, double *, double *, double *, double *, double *);
void rkf45_rvv_4(int, double, double *, double *, double *, double *, double *, double *);
void rkf45_rvv_5(int, double, double *, double *, double *, double *, double *, double *, double *);
void rkf45_rvv_6(int, double, double *, double *, double *, double *);
void rkf45_rvv_err(int, double, double *, double *, double *, double *);

void rk2_rvv_1(int, double, double *, double *, double *);
void rk2_rvv_2(int, double, double *, double *, double *);

void rk4_rvv_1(int, double, double *, double *, double *);
void rk4_rvv_2(int, double, double *, double *, double *);
void rk4_rvv_3(int, double, double *, double *, double *);
void rk4_rvv_4(int, double, double *, double *, double *, double *, double *);

void rkf45_1(int, double, double *, double *, double *);
void rkf45_2(int, double, double *, double *, double *, double *);
void rkf45_3(int, double, double *, double *, double *, double *, double *);
void rkf45_4(int, double, double *, double *, double *, double *, double *, double *);
void rkf45_5(int, double, double *, double *, double *, double *, double *, double *, double *);
void rkf45_6(int, double, double *, double *, double *, double *);
void rkf45_err(int, double, double *, double *, double *, double *, double *);

void rk2_1(int, double, double *, double *, double *);
void rk2_2(int, double, double *, double *, double *);

void rk4_1(int, double, double *, double *, double *);
void rk4_2(int, double, double *, double *, double *);
void rk4_3(int, double, double *, double *, double *);
void rk4_4(int, double, double *, double *, double *, double *, double *);