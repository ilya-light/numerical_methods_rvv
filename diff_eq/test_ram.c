#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main()
{
  system("free -b | grep Mem > ramGraph");

  FILE *f = fopen("cpuBefore", "r");
  FILE *f2 = fopen("cpuAfter", "r");

  unsigned long long oldFirst, oldSecond, oldThird, oldFourth, oldFifth, oldSixth, newFirst, newSecond, newThird, newFourth, newFifth, newSixth;

  fscanf(f, "Mem: %llu %llu %llu %llu %llu %llu", &oldFirst, &oldSecond, &oldThird, &oldFourth, &oldFifth, &oldSixth);
  fscanf(f2, "Mem: %llu %llu %llu %llu %llu %llu", &newFirst, &newSecond, &newThird, &newFourth, &newFifth, &newSixth);
  fclose(f);
  fclose(f2);

  printf("%llu %llu %llu %llu %llu %llu\n", oldFirst, oldSecond, oldThird, oldFourth, oldFifth, oldSixth);
  printf("%llu %llu %llu %llu %llu %llu\n", newFirst, newSecond, newThird, newFourth, newFifth, newSixth);
  printf("%llu\n", newSecond - oldSecond);
}