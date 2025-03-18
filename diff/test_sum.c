
#include <riscv-vector.h>
#include <stdlib.h>
#include <stdio.h>

int main()
{
    double * array = (double*)malloc(sizeof(double) * 20);
    for(int i =0; i<20;i++)
    {
        array[i] = i;
    }
    double error_d = 0;
    double temp_n = 20;
    size_t vlmax = vsetvlmax_e64m4();
    vfloat64m4_t sum_k = vfmv_v_f_f64m4(0, vlmax);
    vfloat64m1_t vec_sum = vfmv_v_f_f64m1(0.0, 1);
    for(size_t vl; temp_n > 0; temp_n -= vl, array+=vl)
    {
        printf("iter\n");
        vl = vsetvl_e64m4(temp_n);
        sum_k = vle_v_f64m4(array, vl);
        vec_sum = vfredsum_vs_f64m4_f64m1(vec_sum, sum_k, vec_sum, vl);
    }
    error_d = vfmv_f_s_f64m1_f64(vec_sum,1);
    printf("%f", error_d);
}