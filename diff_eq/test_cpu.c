#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "rk.h"

void read_cpu_usage(unsigned long long *total, unsigned long long *idle)
{
  FILE *proc_stat = fopen("/proc/stat", "r");
  if (proc_stat == NULL)
  {
    perror("Failed to open /proc/stat");
    exit(EXIT_FAILURE);
  }
  char line[256];
  if (fgets(line, sizeof(line), proc_stat))
  {
    char cpu_label[5];
    unsigned long long user, nice, system, irq, softirq, steal, guest, guest_nice;
    sscanf(line, "%s %llu %llu %llu %llu %llu %llu %llu %llu %llu", cpu_label, &user, &nice, &system, idle, &irq, &softirq, &steal, &guest, &guest_nice);
    *total = user + nice + system + irq + softirq + steal;
  }
  fclose(proc_stat);
}

double calculate_cpu_usage(unsigned long long prev_total, unsigned long long prev_idle, unsigned long long curr_total, unsigned long long curr_idle)
{
  long long total_diff = abs(curr_total - prev_total);
  long long idle_diff = abs(curr_idle - prev_idle);
  if (total_diff == 0) return 0.0;
  return 100.0 * (total_diff - idle_diff) / total_diff;
}





int N = 0;
const int TEST_COUNT = 10;

void f(double t, double y[], double dydt[]) {
    double b1 = 0.02; 
    double c1 = 0.01; 
    dydt[0] = t;
    for (int i = 1; i < N; i++)
    {
        dydt[i] = t * y[i-1];
    }
}


int main() {
    
    FILE *file = fopen("results_cpu_scalar.txt", "w");

    for (int i = 2; i < 1000; i*=2)
    {
        N = i;
        
        double t0 = 0.0;
        double * y0 = (double*)malloc(sizeof(double) * N);
        double * y = (double*)malloc(sizeof(double) * N);
        double t_end = 3.0;
        double h = 0.1;
        for (int i = 0; i < N; i++)
        {
            y0[i] = i;
        }


        unsigned long long prev_total = 0, prev_idle = 0;
        unsigned long long curr_total = 0, curr_idle = 0;

        read_cpu_usage(&prev_total, &prev_idle);
        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            rkf45(&f, t0, y0, y, t_end, h, N);    
        }

        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            rk2(&f, t0, y0, y, t_end, h, N);
        }
        
        for (size_t i = 0; i < TEST_COUNT; i++)
        {
            rk4(&f, t0, y0, y, t_end, h, N);
        }
        read_cpu_usage(&curr_total, &curr_idle);
        printf("CPU Usage: %.2f%%\n", calculate_cpu_usage(prev_total, prev_idle, curr_total, curr_idle));
        
        free(y0);
        free(y);
    }

    fclose(file);
    return 0;
}
