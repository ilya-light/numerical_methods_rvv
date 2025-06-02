#include "rk.h"

#define N 4
#define H 0.1
#define ACCURAСY 1.0E-14

double * y_scalar;
double * y_rvv;
double * k1;
double * k2;
double * k3;
double * k4;
double * k5;
double * temp_y_scalar;
double * temp_y_rvv;

double * y_rvv_p;
double * k1_p;
double * k2_p;
double * k3_p;
double * k4_p;
double * k5_p;
double * temp_y_rvv_p;

int arrays_equal(double * first, double * second)
{
    for (int i = 0; i < N; i++)
    {
        double diff = fabs(first[i] - second[i]);
        if(diff > ACCURAСY)
        {
            printf("diff %e\n", diff);
            return 0;
        }
    }
    return 1;
}

int test_rkf45_1()
{
    rkf45_1(N, H, y_scalar, k1, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    temp_y_rvv_p = temp_y_rvv;
    rkf45_rvv_1(N, H, y_rvv_p, k1_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rkf45_2()
{
    rkf45_2(N, H, y_scalar, k1, k2, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k2_p = k2;
    temp_y_rvv_p = temp_y_rvv;
    rkf45_rvv_2(N, H, y_rvv_p, k1_p, k2_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rkf45_3()
{
    rkf45_3(N, H, y_scalar, k1, k2, k3, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k2_p = k2;
    k3_p = k3;
    temp_y_rvv_p = temp_y_rvv;
    rkf45_rvv_3(N, H, y_rvv_p, k1_p, k2_p, k3_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rkf45_4()
{
    rkf45_4(N, H, y_scalar, k1, k2, k3, k4, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k2_p = k2;
    k3_p = k3;
    k4_p = k4;
    temp_y_rvv_p = temp_y_rvv;
    rkf45_rvv_4(N, H, y_rvv_p, k1_p, k2_p, k3_p, k4_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rkf45_5()
{
    rkf45_5(N, H, y_scalar, k1, k2, k3, k4, k5, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k2_p = k2;
    k3_p = k3;
    k4_p = k4;
    k5_p = k5;
    temp_y_rvv_p = temp_y_rvv;
    rkf45_rvv_5(N, H, y_rvv_p, k1_p, k2_p, k3_p, k4_p, k5_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rkf45_6()
{
    rkf45_6(N, H, y_scalar, k1, k4, k5);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k4_p = k4;
    k5_p = k5;
    rkf45_rvv_6(N, H, y_rvv_p, k1_p, k4_p, k5_p);
    return(arrays_equal(y_rvv, y_scalar));
}

int test_rk2_1()
{
    rk2_1(N, H, y_scalar, k1, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    temp_y_rvv_p = temp_y_rvv;
    rk2_rvv_1(N, H, y_rvv_p, k1_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rk2_2()
{
    rk2_2(N, H, y_scalar, k1, k2);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k2_p = k2;
    rk2_rvv_2(N, H, y_rvv_p, k1_p, k2_p);
    return(arrays_equal(y_rvv, y_scalar));
}

int test_rk4_1()
{
    rk4_1(N, H, y_scalar, k1, temp_y_scalar);

    y_rvv_p = y_rvv;
    k1_p = k1;
    temp_y_rvv_p = temp_y_rvv;
    rk4_rvv_1(N, H, y_rvv_p, k1_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rk4_2()
{
    rk4_2(N, H, y_scalar, k2, temp_y_scalar);

    y_rvv_p = y_rvv;
    k2_p = k2;
    temp_y_rvv_p = temp_y_rvv;
    rk4_rvv_2(N, H, y_rvv_p, k2_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rk4_3()
{
    rk4_3(N, H, y_scalar, k3, temp_y_scalar);

    y_rvv_p = y_rvv;
    k3_p = k3;
    temp_y_rvv_p = temp_y_rvv;
    rk4_rvv_3(N, H, y_rvv_p, k3_p, temp_y_rvv_p);
    return(arrays_equal(temp_y_rvv, temp_y_scalar));
}

int test_rk4_4()
{
    rk4_4(N, H, y_scalar, k1, k2, k3, k4);

    y_rvv_p = y_rvv;
    k1_p = k1;
    k2_p = k2;
    k3_p = k3;
    k4_p = k4;
    rk4_rvv_4(N, H, y_rvv_p, k1_p, k2_p, k3_p, k4_p);

    return(arrays_equal(y_rvv, y_scalar));
}

int test_any(int (*tf)(), char * name)
{
    int result = tf();
    printf("Function %s is tested\nResult:", name);
    if(result == 1)
    {
        printf("success\n");
    }
    else
    {
        printf("fail\n");
    }
}



void init()
{
    y_scalar = (double*)malloc(sizeof(double) * N);
    y_rvv = (double*)malloc(sizeof(double) * N);
    k1 = (double*)malloc(sizeof(double) * N);
    k2 = (double*)malloc(sizeof(double) * N);
    k3 = (double*)malloc(sizeof(double) * N);
    k4 = (double*)malloc(sizeof(double) * N);
    k5 = (double*)malloc(sizeof(double) * N);
    temp_y_scalar = (double*)malloc(sizeof(double) * N);
    temp_y_rvv = (double*)malloc(sizeof(double) * N);
    
    for (int i = 0; i < N; i++)
    {
        k1[i] = i * 2 + 1;
        k2[i] = i * 3 + 1;
        k3[i] = i * 4 + 1;
        k4[i] = i * 5 + 1;
        k5[i] = i * 6 + 1;
        y_scalar[i] = i * 7 + 1;
        y_rvv[i] = i * 7 + 1;
    }
}

void free_arrays()
{
    free(y_rvv);
    free(y_scalar);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(temp_y_rvv);
    free(temp_y_scalar);
}

int main()
{
    init();

    test_any(test_rkf45_1, "rkf45_1");
    test_any(test_rkf45_2, "rkf45_2");
    test_any(test_rkf45_3, "rkf45_3");
    test_any(test_rkf45_4, "rkf45_3");
    test_any(test_rkf45_5, "rkf45_3");
    test_any(test_rkf45_6, "rkf45_3");
    test_any(test_rk4_1, "rk4_1");
    test_any(test_rk4_2, "rk4_2");
    test_any(test_rk4_3, "rk4_3");
    test_any(test_rk4_4, "rk4_4");
    test_any(test_rk2_1, "rk2_1");
    test_any(test_rk2_2, "rk2_2");

    free_arrays();
}

