#include "utilities.h"

void metropolis(int M, int N, int step, double T, double *en, double *mag,
                double *dE, double *spin) {
  /* given configurations matrix spin, updates state of spin #step according to
   * the energy deltas */
  int i;
  double rr = 0.0;
#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for (i = 0; i < M; ++i) {
    rr = (double)rand() / RAND_MAX;
    if (rr < exp(-dE[i] / T)) { // oss: this dE is E_final - E_initial, so we
                                // expect configuation to be interesting if dE<0
      spin[i * N + step] = -spin[i * N + step];
      // update energy
      en[i] += dE[i];
      magnetization(spin, mag, N, i);
      // mag[i]+=2.0/N*spin[i*N+step]; // nb: there might be a - sign here
    }
  }
}

void transpose(double *start, double *end, int N, int M) {
/* transpose start (N*M) in end (M*N) */
#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      end[j * N + i] = start[i * M + j];
    }
  }
}

void average(double *x, int N, double *out) {
  /* returns average value of array x in variable out */
  *out = 0.0;
#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < N; ++i)
    *out += x[i];
  *out /= 1.0 * N;
}

void magnetization(double *spin, double *mag, int N, int idx) {
  double tmp = 0.0;
  for (int j = 0; j < N; j++)
    tmp = tmp + spin[idx * N + j];
  tmp = tmp / N;
  mag[idx] = tmp * tmp;
}

double cclock() {
  struct timeval tmp;
  double sec;
  gettimeofday(&tmp, (struct timezone *)0);
  sec = tmp.tv_sec + ((double)tmp.tv_usec) / 1000000.0;
  return sec;
}

void stder(double T, double *data, int N, FILE *file) {
  double sum = 0.0, mean, stder = 0.0;
  double sum1 = 0.0, mean1, stder1 = 0.0;
  int i;
  // #ifdef _USE_OMP
  // #pragma omp parallel for
  // #endif
  for (i = 0; i < N; ++i) {
    if (i % 2 == 0) {
      sum = sum + data[i];
    } else {
      sum1 = sum1 + data[i];
    }
  }
  mean = sum * 2 / N;
  mean1 = sum1 * 2 / N;
  // #ifdef _USE_OMP
  // #pragma omp parallel for
  // #endif
  for (i = 0; i < N; ++i) {
    if (i % 2 == 0) {
      stder += pow(data[i] - mean, 2);
    } else {
      stder1 += pow(data[i] - mean1, 2);
    }
  }
  int ii = 0;
  // for (int kk = 0; kk < N; kk += 2) {
  //   printf("%d %f %f\n", ii++, data[kk], data[kk + 1]);
  // }
  stder = sqrt(stder / (2 * N));
  stder1 = sqrt(stder1 / (2 * N));

  printf("%f %f %f %f %f\n", T, mean, stder, mean1, stder1);
}
