#include "math.h"
#include "stdio.h"
#include "unistd.h"

int lcg(int m, int a, int c, int x0);

int main(int argc, char *argv[]) {
  unsigned int x = 16807;      /* 7^5 */
  unsigned int m = 2147483647; /* 2^31 -1 */
  lcg(m, x, 3, 1);
  // lcg(m, x, 0, 1); // mlcg random number generator
  return 0;
}

int lcg(int m, int a, int c, int x0) {
  int x;
  while (1) {
    x = a * x0 + c;
    x0 = x % m;
    write(1, &x0, 4);
  }
  return 0;
}
