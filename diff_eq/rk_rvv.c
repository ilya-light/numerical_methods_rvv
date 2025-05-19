#include <riscv-vector.h>
#include "rk.h"

void copy_vector_rvv(double * src, double * dest, int n)
{
    for(size_t vl; n > 0; n -= vl, src+=vl, dest+=vl)
    {
        vl = vsetvl_e64m4(n);
        vfloat64m4_t vec_src = vle_v_f64m4(src,vl);
        vse_v_f64m4(dest, vec_src, vl);
    }
}


size_t vlmax;
vfloat64m4_t k1_vec;
vfloat64m4_t k2_vec;
vfloat64m4_t k3_vec;
vfloat64m4_t k4_vec;
vfloat64m4_t k5_vec;
vfloat64m4_t k6_vec; 
vfloat64m4_t y_vec;
vfloat64m4_t sum_k;
vfloat64m1_t vec_sum;

void init_rvv()
{
    vlmax = vsetvlmax_e64m4();
    k1_vec = vfmv_v_f_f64m4(0, vlmax);
    k2_vec = vfmv_v_f_f64m4(0, vlmax);
    k3_vec = vfmv_v_f_f64m4(0, vlmax);
    k4_vec = vfmv_v_f_f64m4(0, vlmax);
    k5_vec = vfmv_v_f_f64m4(0, vlmax);
    k6_vec = vfmv_v_f_f64m4(0, vlmax);
    y_vec = vfmv_v_f_f64m4(0, vlmax);
    sum_k = vfmv_v_f_f64m4(0, vlmax);
    vec_sum = vfmv_v_f_f64m1(0.0, 1);
}

void rkf45_rvv(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n)
{
    init_rvv();
    copy_vector_rvv(y0, y, n);


    double *k1 = (double*)malloc(n * sizeof(double));
    double *k2 = (double*)malloc(n * sizeof(double));
    double *k3 = (double*)malloc(n * sizeof(double));
    double *k4 = (double*)malloc(n * sizeof(double));
    double *k5 = (double*)malloc(n * sizeof(double));
    double *k6 = (double*)malloc(n * sizeof(double));
    double *y_temp = (double*)malloc(n * sizeof(double));
  
    double * k1_p, *k2_p, *k3_p, *k4_p, *k5_p, *k6_p, *y_p, *temp_y_p, *err_p;
    int temp_n;
    double t = t0;

    while (t < t_end)
    {
        f(t, y, k1);

        k1_p = k1;
        y_p = y;
        temp_y_p = y_temp;
        temp_n = n;

        rkf45_rvv_1(temp_n, h, y_p, k1_p, temp_y_p);

        f(t + h / 4, y_temp, k2);
        
        k1_p = k1;
        k2_p = k2;
        y_p = y;
        temp_y_p = y_temp;
        temp_n = n;
        
        rkf45_rvv_2(temp_n, h, y_p, k1_p, k2_p, temp_y_p);

        f(t + 3 * h / 8, y_temp, k3);
        
        k1_p = k1;
        k2_p = k2;
        k3_p = k3;
        y_p = y;
        temp_y_p = y_temp;
        temp_n = n;

        rkf45_rvv_3(temp_n, h, y_p, k1_p, k2_p, k3_p, temp_y_p);

        f(t + 12 * h / 13, y_temp, k4);
        
        k1_p = k1;
        k2_p = k2;
        k3_p = k3;
        k4_p = k4;
        y_p = y;
        temp_y_p = y_temp;
        temp_n = n;
        
        rkf45_rvv_4(temp_n, h, y_p, k1_p, k2_p, k3_p, k4_p, temp_y_p);

        f(t + h, y_temp, k5);
        
        k1_p = k1;
        k2_p = k2;
        k3_p = k3;
        k4_p = k4;
        k5_p = k5;
        y_p = y;
        temp_y_p = y_temp;
        temp_n = n;
        
        rkf45_rvv_5(temp_n, h, y_p, k1_p, k2_p, k3_p, k4_p, k5_p, temp_y_p);

        f(t + h / 2, y_temp, k6);

        k1_p = k1;
        k4_p = k4;
        k5_p = k5;
        y_p = y;
        temp_n = n;
        
        rkf45_rvv_6(temp_n, h, y_p, k1_p, k4_p, k5_p);        

        k1_p = k1;
        k2_p = k2;
        k3_p = k3;
        k4_p = k4;
        temp_n = n;

        double error_d = 0;
        rkf45_rvv_err(temp_n, h, k1_p, k2_p,k3_p, k4_p);
        error_d = vfmv_f_s_f64m1_f64(vec_sum,1);
        error_d = sqrt(error_d);

        if (error_d > 1e-6)
        {
            h *= 0.5;
        }
        else if (error_d < 1e-8)
        {
            h *= 1.5;
        }
        t += h;
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(y_temp);
}

void rkf45_rvv_1(int temp_n, double h,  double * y_p, double * k1_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, h/4, k1_vec, vl);
            
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rkf45_rvv_2(int temp_n, double h, double * y_p, double * k1_p, double * k2_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k2_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        k2_vec = vle_v_f64m4(k2_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, 3*h/32, k1_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, 9*h/32, k2_vec, vl);
            
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rkf45_rvv_3(int temp_n, double h, double * y_p, double * k1_p, double * k2_p, double * k3_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k2_p+=vl, k3_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        k2_vec = vle_v_f64m4(k2_p, vl);
        k3_vec = vle_v_f64m4(k3_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, 1932*h/2197, k1_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, -7200*h/2197, k2_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, 7296*h/2197, k3_vec, vl);
            
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rkf45_rvv_4(int temp_n, double h, double * y_p, double * k1_p, double * k2_p, double * k3_p, double * k4_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k2_p+=vl, k3_p+=vl, k4_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        k2_vec = vle_v_f64m4(k2_p, vl);
        k3_vec = vle_v_f64m4(k3_p, vl);
        k4_vec = vle_v_f64m4(k4_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, 439*h/216, k1_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, -8*h/27, k2_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, 3680*h/513, k3_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, -845*h/4104, k4_vec, vl);
            
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rkf45_rvv_5(int temp_n, double h, double * y_p, double * k1_p, double * k2_p, double * k3_p, double * k4_p, double * k5_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k2_p+=vl, k3_p+=vl, k4_p+=vl, k5_p+=vl, temp_y_p+=vl)
{
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        k2_vec = vle_v_f64m4(k2_p, vl);
        k3_vec = vle_v_f64m4(k3_p, vl);
        k4_vec = vle_v_f64m4(k4_p, vl);
        k5_vec = vle_v_f64m4(k5_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, -8*h/27, k1_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, 2*h, k2_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, -3544*h/2565, k3_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, 1859*h/4104, k4_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, -11*h/40, k5_vec, vl);
            
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rkf45_rvv_6(int temp_n, double h, double * y_p, double * k1_p, double * k4_p, double * k5_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k4_p+=vl, k5_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        k4_vec = vle_v_f64m4(k4_p, vl);
        k5_vec = vle_v_f64m4(k5_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        sum_k = vfadd_vv_f64m4(k1_vec, k5_vec, vl);
        sum_k = vfmacc_vf_f64m4(sum_k, 4, k4_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, h/360, sum_k, vl);
            
        vse_v_f64m4(y_p, y_vec, vl);
    }
}

void rkf45_rvv_err(int temp_n, double h, double * k1_p, double * k2_p, double * k3_p, double * k4_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, k1_p+=vl, k2_p+=vl, k3_p+=vl, k4_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        k2_vec = vle_v_f64m4(k2_p, vl);
        k3_vec = vle_v_f64m4(k3_p, vl);
        k4_vec = vle_v_f64m4(k4_p, vl);
            
        sum_k = vfsub_vv_f64m4(k1_vec, k4_vec, vl);
        sum_k = vfmacc_vf_f64m4(sum_k, -2, k2_vec, vl);
        sum_k = vfmacc_vf_f64m4(sum_k, 2, k3_vec, vl);
        sum_k = vfmul_vf_f64m4(sum_k, h/360, vl);
        sum_k = vfmul_vv_f64m4(sum_k, sum_k, vl);

        vec_sum = vfredsum_vs_f64m4_f64m1(vec_sum, sum_k, vec_sum, vl);
    }
}

void rk2_rvv(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n)
{    
    init_rvv();
    copy_vector_rvv(y0, y, n);
    
    double * k1 = (double*)malloc(n * sizeof(double));
    double * k2 = (double*)malloc(n * sizeof(double));
    double * y_temp = (double*)malloc(n * sizeof(double));

    double * k1_p, *k2_p, *y_p, *temp_y_p;
    double t = t0;
    double hs = h / 2;
    double temp_n;

    while (t < t_end)
    {
        f(t, y, k1);

        temp_n = n;
        k1_p = k1;
        y_p = y;
        temp_y_p = y_temp;

        rk2_rvv_1(temp_n, hs, y_p, k1_p, temp_y_p);

        f(t + h, y_temp, k2);

        temp_n = n;
        k1_p = k1;
        k2_p = k2;
        y_p = y;
        
        rk2_rvv_2(temp_n, hs, y_p, k1_p, k2_p);

        t += h;

    }
    free(k1);
    free(k2);
    free(y_temp);
}

void rk2_rvv_1(int temp_n, double hs, double * y_p, double * k1_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);

        k1_vec = vle_v_f64m4(k1_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, hs, k1_vec, vl);
            
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rk2_rvv_2(int temp_n, double hs, double * y_p, double * k1_p, double * k2_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k2_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);
            
        k1_vec = vle_v_f64m4(k1_p, vl); 
        k2_vec = vle_v_f64m4(k2_p, vl);
        y_vec = vle_v_f64m4(y_p, vl);
            
        sum_k = vfadd_vv_f64m4(k1_vec, k2_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, hs, sum_k, vl);

        vse_v_f64m4(y_p, y_vec, vl);
    }
}

void rk4_rvv(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n)
{    
    init_rvv();
    copy_vector_rvv(y0, y, n);

    double * k1 = (double*)malloc(n * sizeof(double));
    double * k2 = (double*)malloc(n * sizeof(double));
    double * k3 = (double*)malloc(n * sizeof(double));
    double * k4 = (double*)malloc(n * sizeof(double));
    double * y_temp = (double*)malloc(n * sizeof(double));
    
    double * k1_p, *k2_p, *k3_p, *k4_p, *y_p, *temp_y_p;
    double t = t0;
    double temp_n, hs;
    
    while (t < t_end)
    {
        hs = h/2;

        f(t, y, k1);
        
        temp_n = n;
        y_p = y; 
        k1_p = k1;
        temp_y_p = y_temp;
        
        rk4_rvv_1(temp_n, hs, y_p, k1_p, temp_y_p);

        f(t + hs, y_temp, k2);
        
        temp_n = n;
        y_p = y; 
        k2_p = k2;
        temp_y_p = y_temp;

        rk4_rvv_2(temp_n, hs, y_p, k2_p, temp_y_p);
        
        f(t + hs, y_temp, k3);

        temp_n = n;
        y_p = y; 
        k3_p = k3;
        temp_y_p = y_temp;

        rk4_rvv_3(temp_n, h, y_p, k3_p, temp_y_p);

        f(t + h, y_temp, k4);

        hs /= 3;
        temp_n = n;
        y_p = y; 
        k1_p = k1;
        k2_p = k2;
        k3_p = k3;
        k4_p = k4;

        rk4_rvv_4(temp_n, hs, y_p, k1_p, k2_p, k3_p, k4_p);

        t += h;

    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(y_temp);
}

void rk4_rvv_1(int temp_n, double hs, double * y_p, double * k1_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);
        k1_vec = vle_v_f64m4(k1_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, hs, k1_vec, vl);
        
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rk4_rvv_2(int temp_n, double hs, double * y_p, double * k2_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k2_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);
        k2_vec = vle_v_f64m4(k2_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, hs, k2_vec, vl);
        
        vse_v_f64m4(temp_y_p, y_vec, vl);
        }
}

void rk4_rvv_3(int temp_n, double h, double * y_p, double * k3_p, double * temp_y_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k3_p+=vl, temp_y_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);
        k3_vec = vle_v_f64m4(k3_p, vl);
        y_vec = vle_v_f64m4(y_p,vl);
            
        y_vec = vfmacc_vf_f64m4(y_vec, h, k3_vec, vl);
        
        vse_v_f64m4(temp_y_p, y_vec, vl);
    }
}

void rk4_rvv_4(int temp_n, double hs, double * y_p, double * k1_p, double * k2_p, double * k3_p, double * k4_p)
{
    for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k1_p+=vl, k2_p+=vl, k3_p+=vl, k4_p+=vl)
    {
        vl = vsetvl_e64m4(temp_n);
            
        k1_vec = vle_v_f64m4(k1_p, vl); 
        k2_vec = vle_v_f64m4(k2_p, vl);
        k3_vec = vle_v_f64m4(k3_p, vl);
        k4_vec = vle_v_f64m4(k4_p, vl);
        y_vec = vle_v_f64m4(y_p, vl);
            
        sum_k = vfadd_vv_f64m4(k1_vec, k4_vec, vl);
        sum_k = vfmacc_vf_f64m4(sum_k, 2, k2_vec, vl);
        sum_k = vfmacc_vf_f64m4(sum_k, 2, k3_vec, vl);
        y_vec = vfmacc_vf_f64m4(y_vec, hs, sum_k, vl);

        vse_v_f64m4(y_p, y_vec, vl);
    }
}
