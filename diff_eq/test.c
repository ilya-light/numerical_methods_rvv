#include <time.h>
#include <math.h>
#include "rk.h"

int N = 0;
const int TEST_COUNT = 30;

void f(double t, double y[], double dydt[]) {
    double b1 = 0.02; 
    double c1 = 0.01; 
    dydt[0] = t;
    for (int i = 1; i < N; i++)
    {
        dydt[i] = t * y[i-1];
    }
}

int main()
{
    FILE *file = fopen("results_vector.txt", "w");

    for (int i = 2; i < 1000; i*=2)
    {
        N = i;
        
        double t0 = 0.0;
        double * y0 = (double*)malloc(sizeof(double) * N);
        double * y = (double*)malloc(sizeof(double) * N);
        double t_end = 5.0;
        double h = 0.1;
        for (int i = 0; i < N; i++)
        {
            y0[i] = i;
        }

        double sum = 0;

        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            clock_t start = clock();
            rkf45(&f, t0, y0, y, t_end, h, N);
            clock_t end = clock();
            sum += ((double)(end - start))/CLOCKS_PER_SEC;
        }
        fprintf(file, "rkf45,%d,%.6f\n", N, sum / TEST_COUNT);
        
        sum = 0;
        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            clock_t start = clock();
            rk2(&f, t0, y0, y, t_end, h, N);
            clock_t end = clock();
            sum += ((double)(end - start))/CLOCKS_PER_SEC;
        }
        fprintf(file, "rk2,%d,%.6f\n", N, sum / TEST_COUNT);
        
        sum = 0;
        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            clock_t start = clock();
            rk4(&f, t0, y0, y, t_end, h, N);
            clock_t end = clock();
            sum += ((double)(end - start))/CLOCKS_PER_SEC;
        }
        fprintf(file, "rk4,%d,%.6f\n", N, sum / TEST_COUNT);
        
        free(y0);
        free(y);
    }

    fclose(file);
    return 0;
}