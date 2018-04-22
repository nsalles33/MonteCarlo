#include "math.h"
#include "stdio.h"
#include "unistd.h"

int bit_gen(int m);

int main(int argc, char *argv[]) {
  // unsigned int x = 16807;      /* 7^5 */
  // unsigned int m = 2147483647; /* 2^31 -1 */
  bit_gen(15);
  // lcg(m, x, 0, 1); // mlcg random number generator
  return 0;
}

int bit_gen(int a) {
  unsigned int x;
  while (1) {
    x = (((a & 8) >> 3) ^ (a & 1)) << 5;
    a = a >> 1;
    a = a | x;
    // printf("%i\n", a);
    write(1, &a, 4);
  }
  return 0;
}
