#include "rk4.h"


// Реализация метода Рунге-Кутты 4-го порядка
void runge_kutta(void (*f)(double, double*, double*), double t0, double * y0, double t_end, double h, int n) {
    double t = t0;
    double * y = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        y[i] = y0[i];
    }

    while (t < t_end) {
        double * k1 = (double*)malloc(n * sizeof(double));
        double * k2 = (double*)malloc(n * sizeof(double));
        double * k3 = (double*)malloc(n * sizeof(double));
        double * k4 = (double*)malloc(n * sizeof(double));
        double * y_temp = (double*)malloc(n * sizeof(double));

        f(t, y, k1); // k1 = f(t, y)
        for (int i = 0; i < n; i++) {
            y_temp[i] = y[i] + 0.5 * h * k1[i];
        }
        f(t + 0.5 * h, y_temp, k2); // k2 = f(t + h/2, y + h/2 * k1)
        for (int i = 0; i < n; i++) {
            y_temp[i] = y[i] + 0.5 * h * k2[i];
        }
        f(t + 0.5 * h, y_temp, k3); // k3 = f(t + h/2, y + h/2 * k2)
        for (int i = 0; i < n; i++) {
            y_temp[i] = y[i] + h * k3[i];
        }
        f(t + h, y_temp, k4); // k4 = f(t + h, y + h * k3)

        for (int i = 0; i < n; i++) {
            y[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
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

