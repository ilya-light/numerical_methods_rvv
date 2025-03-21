#include <time.h>
#include <math.h>
#include "rk.h"

#define FILENAME_SIZE 30

int N = 0;
const int TEST_COUNT = 100;

struct method_params_t
{
    void (*f)(double, double *, double *);
    double t0;
    double * y0;
    double * y;
    double t_end;
    double h;
};

void f(double t, double y[], double dydt[]) {
    dydt[0] = t;
    for (int i = 1; i < N; i++)
    {
        dydt[i] = t * y[i-1];
    }
}

double test_method(void (*method)(void (*)(double, double *, double *), double, double *, double *, double, double, int),
                    struct method_params_t params)
{
    double sum = 0;

    for (size_t i = 0; i < TEST_COUNT; i++)
    {
        clock_t start = clock();
        method(params.f, params.t0, params.y0, params.y, params.t_end, params.h, N);
        clock_t end = clock();
        sum += ((double)(end - start))/CLOCKS_PER_SEC;
    }
    return sum / TEST_COUNT;
}

int main(int argc, char *argv[])
{
    char filename[FILENAME_SIZE]; 
    if(argc > 1 && argv[1] == "-v")
    {
        strncpy(filename, "results_vector.txt", FILENAME_SIZE);
    }
    else
    {
        strncpy(filename, "results_scalar.txt", FILENAME_SIZE);
    }
    FILE *file = fopen(filename, "w");

    for (int i = 2; i <= 512; i*=2)
    {
        N = i;
        
        struct method_params_t params;
        params.t0 = 0.0;
        params.y0 = (double*)malloc(sizeof(double) * N);
        params.y = (double*)malloc(sizeof(double) * N);
        params.t_end = 5.0;
        params.h = 0.1;
        for (int i = 0; i < N; i++)
        {
            params.y0[i] = i;
        }

        fprintf(file, "rkf45,%d,%.6f\n", N, test_method(&rkf45, params));
        fprintf(file, "rk2,%d,%.6f\n", N, test_method(&rkf45, params));
        fprintf(file, "rk4,%d,%.6f\n", N, test_method(&rkf45, params));
        
        free(params.y0);
        free(params.y);
    }

    fclose(file);
    return 0;
}