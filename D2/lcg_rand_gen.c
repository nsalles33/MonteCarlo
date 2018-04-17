#include "math.h"
#include "stdio.h"
#include "unistd.h"

int lcg(int m, int a, int c, int x0, int n);

int main(int argc, char *argv[]) {
  lcg(83, 7, 3, 1, 30);
  // lcg(83, 7, 0, 1, 30); // mlcg random number generator
  return 0;
}

int lcg(int m, int a, int c, int x0, int n) {
  int x;
  for (int i = 0; i < n; i++) {
    x = a * x0 + c;
    x0 = x % m;
    printf("%i\n", x0);
  }
  return 0;
}
