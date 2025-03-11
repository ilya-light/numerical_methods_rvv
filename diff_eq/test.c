#include "rk4.h"

#include <math.h>

const int n = 50000;

// Определение функций роста, зависящих от времени
double a1(double t) {
    return 0.1 + 0.05 * sin(0.1 * t); // Пример функции роста для первого вида
}

// Определение функции, представляющей систему ОДУ
void f(double t, double y[], double dydt[]) {
    double b1 = 0.02; // Коэффициент взаимодействия для первого вида
    double c1 = 0.01; // Влияние третьего вида на первый

    for (int i = 0; i < n; i++)
    {
        dydt[i] = a1(t) * y[0] - b1 * y[i] * y[1] + c1 * y[2];
    }
}

int main() {
    double t0 = 0.0; // Начальное время
    double * y0 = (double*)malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
    {
        y0[i] = i;
    }
    
    double t_end = 50.0; // Конечное время
    double h = 0.1; // Шаг интегрирования

    runge_kutta(&f, t0, y0, t_end, h, n); // Вызов функции Рунге-Кутты
    return 0;
}