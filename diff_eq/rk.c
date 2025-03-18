#include <time.h>
#include "rk.h"

void copy_vector(double * src, double * dest, int n)
{
    for (int i = 0; i < n; i++)
    {
        dest[i] = src[i];
    }
}

void rkf45(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n)
{
    copy_vector(y0, y, n);

    double *k1 = (double*)malloc(n * sizeof(double));
    double *k2 = (double*)malloc(n * sizeof(double));
    double *k3 = (double*)malloc(n * sizeof(double));
    double *k4 = (double*)malloc(n * sizeof(double));
    double *k5 = (double*)malloc(n * sizeof(double));
    double *k6 = (double*)malloc(n * sizeof(double));
    double *y_temp = (double*)malloc(n * sizeof(double));

    double t = t0;

    while (t < t_end)
    {
        f(t, y, k1);
        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + (h / 4) * k1[i];
        }
        f(t + h / 4, y_temp, k2);
        
        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + (3 * h / 32) * k1[i] + (9 * h / 32) * k2[i];
        }
        f(t + 3 * h / 8, y_temp, k3);
        
        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + (1932 * h / 2197) * k1[i] - (7200 * h / 2197) * k2[i] + (7296 * h / 2197) * k3[i];
        }
        f(t + 12 * h / 13, y_temp, k4);
        
        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + (439 * h / 216) * k1[i] - (8 * h / 27) * k2[i] + (3680 * h / 513) * k3[i] - (845 * h / 4104) * k4[i];
        }
        f(t + h, y_temp, k5);
        
        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] - (8 * h / 27) * k1[i] + 2 * h * k2[i] - (3544 * h / 2565) * k3[i] + (1859 * h / 4104) * k4[i] - (11 * h / 40) * k5[i];
        }
        f(t + h / 2, y_temp, k6);

        for (int i = 0; i < n; i++)
        {
            y[i] += (h / 360) * (k1[i] + 4 * k4[i] + k5[i]);
        }

        double error = 0.0;
        for (int i = 0; i < n; i++)
        {
            double temp_error = (h / 360) * (k1[i] - 2 * k2[i] + 2 * k3[i] - k4[i]);
            error += temp_error * temp_error;
        }
        error = sqrt(error);

        if (error > 1e-6)
        {
            h *= 0.5;
        }
        else if (error < 1e-8)
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

void rk2(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n)
{    
    copy_vector(y0, y, n);
    
    double * k1 = (double*)malloc(n * sizeof(double));
    double * k2 = (double*)malloc(n * sizeof(double));
    double * y_temp = (double*)malloc(n * sizeof(double));

    double t = t0;
    double hs = h / 2;
    while (t < t_end)
    {
        f(t, y, k1);

        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + hs * k1[i];
        }
        f(t + h, y_temp, k2);
        for (int i = 0; i < n; i++)
        {
            y[i] += hs * (k1[i] + k2[i]);
        }
        t += h;
    }
    free(k1);
    free(k2);
    free(y_temp);
}

void rk4(void (*f)(double, double*, double*), double t0, double * y0, double * y, double t_end, double h, int n)
{    
    copy_vector(y0, y, n);

    double * k1 = (double*)malloc(n * sizeof(double));
    double * k2 = (double*)malloc(n * sizeof(double));
    double * k3 = (double*)malloc(n * sizeof(double));
    double * k4 = (double*)malloc(n * sizeof(double));
    double * y_temp = (double*)malloc(n * sizeof(double));
    
    double t = t0;
    double hs;
    
    while (t < t_end)
    {
        hs = h/2;

        f(t, y, k1);

        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + hs * k1[i];
        }

        f(t + hs, y_temp, k2);

        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + hs * k2[i];
        }

        f(t + hs, y_temp, k3);

        for (int i = 0; i < n; i++)
        {
            y_temp[i] = y[i] + h * k3[i];
        }

        hs /= 3;
        f(t + h, y_temp, k4);

        for (int i = 0; i < n; i++)
        {
            y[i] += hs * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }

        t += h;

    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(y_temp);
}
