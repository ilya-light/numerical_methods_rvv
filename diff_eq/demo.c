#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
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

void test_method(void (*method)(void (*)(double, double *, double *), double, double *, double *, double, double, int),
                    struct method_params_t* params,
                    FILE* file,
                    const char* method_name)
{
    for (size_t i = 0; i < TEST_COUNT; i++)
    {
        clock_t start = clock();
        method(params->f, params->t0, params->y0, params->y, params->t_end, params->h, N);
        clock_t end = clock();
        fprintf(file, "%s,%d,%.6f\n", method_name, N, ((double)(end - start))/CLOCKS_PER_SEC);
    }
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

    for (int i = 512; i <= 512; i*=2)
    {
        N = i;
        
        struct method_params_t* params = (struct method_params_t*)malloc(sizeof(struct method_params_t));
        params->f = &f;
        params->t0 = 0.0;
        params->y0 = (double*)malloc(sizeof(double) * N);
        params->y = (double*)malloc(sizeof(double) * N);
        params->t_end = 3.0;
        params->h = 0.1;
        for (int i = 0; i < N; i++)
        {
            params->y0[i] = i;
        }

        test_method(&rkf45, params, file, "rkf45");
        printf("rkf45 N %d done\n", N);
        test_method(&rk4, params, file, "rk4");
        printf("rkf4 N %d done\n", N);
        test_method(&rk2, params, file, "rk2");
        printf("rkf2 N %d done\n", N);
        
        free(params->y0);
        free(params->y);
        free(params);
    }

    fclose(file);
    return 0;
}
