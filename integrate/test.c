#include <stdio.h>
#include <stdlib.h>
#include "trapetional_integrate.h"

double my_function(double x) {
    return x * x;
}

int main() {
    double low = 0;
    double up = 1;
    size_t calls = 10000000;

    
    double result = trapetional_integrate(&my_function, 0, 1, calls);

    // Выводим результат
    printf("Интеграл: %.6f\n", result);

    return 0;
}