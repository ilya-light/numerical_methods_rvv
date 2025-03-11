#include "rk4.h"
#include <riscv-vector.h>
#include <stdio.h>

// Реализация метода Рунге-Кутты 4-го порядка
void runge_kutta(void (*f)(double, double*, double*), double t0, double * y0, double t_end, double h, int n) {
    double t = t0;
    double hs = h/2;

    double * y = (double*)malloc(n * sizeof(double));
    
    int temp_n = n;
    double * src = y0;
    double * dst = y;
    for(size_t vl; temp_n > 0; temp_n -= vl, src+=vl, dst+=vl)
    {
        vl = vsetvl_e64m4(temp_n);
        vfloat64m4_t vec_src = vle_v_f64m4(src,vl);
        vse_v_f64m4(dst, vec_src, vl);
    }

    while (t < t_end) {
        double * k1 = (double*)malloc(n * sizeof(double));
        double * k2 = (double*)malloc(n * sizeof(double));
        double * k3 = (double*)malloc(n * sizeof(double));
        double * k4 = (double*)malloc(n * sizeof(double));
        double * y_temp = (double*)malloc(n * sizeof(double));


        temp_n = n;
        double * y_p = y; 
        double * k_p = k1;
        dst = y_temp;
        
        f(t, y, k1);
        for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k_p+=vl, dst+=vl)
        {
            vl = vsetvl_e64m4(temp_n);
            vfloat64m4_t k1_vec = vle_v_f64m4(k_p, vl); // загружаем массив k1
            vfloat64m4_t hs_k1_vec = vfmul_vf_f64m4(k1_vec, hs, vl); // делаем hs * k1
            
            vfloat64m4_t y_vec = vle_v_f64m4(y_p,vl); // загружаем массив y
            vfloat64m4_t y_temp_vec = vfadd_vv_f64m4(y_vec, hs_k1_vec, vl); // делаем сложение
        
            vse_v_f64m4(dst, y_temp_vec, vl);
        }
        
        temp_n = n;
        y_p = y; 
        k_p = k2;
        dst = y_temp;

        f(t + 0.5 * h, y_temp, k2); // k2 = f(t + h/2, y + h/2 * k1)
        for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k_p+=vl, dst+=vl)
        {
            vl = vsetvl_e64m4(temp_n);
            vfloat64m4_t k2_vec = vle_v_f64m4(k_p, vl); // загружаем массив k1
            vfloat64m4_t hs_k2_vec = vfmul_vf_f64m4(k2_vec, hs, vl); // делаем hs * k1
            
            vfloat64m4_t y_vec = vle_v_f64m4(y_p,vl); // загружаем массив y
            vfloat64m4_t y_temp_vec = vfadd_vv_f64m4(y_vec, hs_k2_vec, vl); // делаем сложение
        
            vse_v_f64m4(dst, y_temp_vec, vl);
        }

        temp_n = n;
        y_p = y; 
        k_p = k3;
        dst = y_temp;

        f(t + 0.5 * h, y_temp, k3);
        for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k_p+=vl, dst+=vl)
        {
            vl = vsetvl_e64m4(temp_n);
            vfloat64m4_t k3_vec = vle_v_f64m4(k_p, vl); // загружаем массив k1
            vfloat64m4_t hs_k3_vec = vfmul_vf_f64m4(k3_vec, hs, vl); // делаем hs * k1
            
            vfloat64m4_t y_vec = vle_v_f64m4(y_p,vl); // загружаем массив y
            vfloat64m4_t y_temp_vec = vfadd_vv_f64m4(y_vec, hs_k3_vec, vl); // делаем сложение
        
            vse_v_f64m4(dst, y_temp_vec, vl);
        }


        hs /= 3;
        temp_n = n;
        y_p = y; 
        k_p = k1;
        double * k2_p = k2;
        double * k3_p = k3;
        double * k4_p = k4;

        f(t + h, y_temp, k4);
        for(size_t vl; temp_n > 0; temp_n -= vl, y_p+=vl, k_p+=vl, k2_p+=vl, k3_p+=vl, k4_p+=vl)
        {
            vl = vsetvl_e64m4(temp_n);
            
            vfloat64m4_t k1_vec = vle_v_f64m4(k_p, vl); 
            vfloat64m4_t k2_vec = vle_v_f64m4(k2_p, vl);
            vfloat64m4_t k3_vec = vle_v_f64m4(k3_p, vl);
            vfloat64m4_t k4_vec = vle_v_f64m4(k4_p, vl);
             
            k2_vec = vfmul_vf_f64m4(k2_vec, 2, vl);
            k3_vec = vfmul_vf_f64m4(k3_vec, 2, vl);
            
            vfloat64m4_t y_vec = vfadd_vv_f64m4(k1_vec, k2_vec, vl);
            y_vec = vfadd_vv_f64m4(y_vec, k3_vec, vl);
            y_vec = vfadd_vv_f64m4(y_vec, k4_vec, vl);
            y_vec = vfmul_vf_f64m4(y_vec, hs, vl);

            vse_v_f64m4(y_p, y_vec, vl);
        }

        t += h;

        free(k1);
        free(k2);
        free(k3);
        free(k4);
        free(y_temp);
    }

    // Вывод результатов
    //for (int i = 0; i < n; i++) {
    //    printf("y[%d] = %f\n", i, y[i]);
    //}

    free(y);
}

