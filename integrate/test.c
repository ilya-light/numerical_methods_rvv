#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "trapetional_integrate.h"

#define TEST_COUNT 100

double my_function(double x) {
    return x * x;
}

int main() {
    double low = 0;
    double up = 1;
    size_t calls = 10000000;

    double sum = 0;
    double result = 0;
        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            clock_t start = clock();
            result = trapetional_integrate(&my_function, 0, 1, calls);
            clock_t end = clock();
            sum += ((double)(end - start))/CLOCKS_PER_SEC;
        }
        printf("trapetional  %.6f\n", sum / TEST_COUNT);    
    printf("Интеграл: %.6f\n", result);

    return 0;
}