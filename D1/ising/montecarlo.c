#include <utilities.h> //all the functions are there

#define M_PI 3.14159265
// #define N 1000     // Number of spins
#define Nsweep 20  // Number of sweeps
#define walker 100 // Number of walkers
#define gamma 1.25
#define NT 1
#define dT 0.1

int main() {

  int i, j;
  double *J, *Ei, *En, *dE, *Mi, *Mf;
  double en = 0., mag = 0.;
  double t1, t2;

  En = (double *)calloc(walker, sizeof(double));
  Ei = (double *)calloc(walker, sizeof(double));
  dE = (double *)calloc(walker, sizeof(double));
  Mi = (double *)calloc(walker, sizeof(double));
  Mf = (double *)calloc(walker, sizeof(double));

  for (int N = 100; N < 101; N = N + 100) {

    double *spin = (double *)calloc(N * walker, sizeof(double)),
           *spinT = (double *)calloc(N * walker, sizeof(double)),
           *spini = (double *)calloc(N * walker, sizeof(double));

    // initializing the system of spin configuration
    J = (double *)calloc(N * N, sizeof(double));
#ifdef _USE_OMP
#pragma omp parallel for
#endif
    for (i = 0; i < N; ++i) {
      for (j = 0; j < walker; ++j) {
        double r = (double)rand() / ((double)RAND_MAX);
        if (r > 0.5)
          spini[j * N + i] = 1.0;
        else
          spini[j * N + i] = -1.0;
      }
    }

    // Calculating the correlation factor J for the configuration
    double a = M_PI / N;
    double b = N / M_PI;
#ifdef _USE_OMP
#pragma omp parallel for
#endif
    for (i = 0; i < N; ++i) {
      for (j = 0; j < N; ++j) {
        if (j == i) {
          J[i * N + j] = 0;
        } else {
          J[i * N + j] = -pow(b * fabs(sin(a * (i - j))), -gamma);
        }
      }
    }

    // Calculating the intial magnetization of the configuration
    double *tmp = (double *)calloc(N, sizeof(double));
    for (i = 0; i < walker; ++i) {
      cblas_dsymv(CblasRowMajor, CblasUpper, N, 1.0, J, N, spini, 1, 0.0, tmp,
                  1);
      Ei[i] = cblas_ddot(N, tmp, 1, spini, 1);
      mag = 0.0;
      // #ifdef _USE_OOMP
      // #pragma omp parallel for
      // #endif
      for (j = 0; j < N; ++j) {
        mag = mag + spini[i * N + j];
      }
      mag = mag / N;
      Mi[i] = mag * mag;
    }

    double T = 10., *buffer;
    int buff = 0;
    // double data[2 * Nsweep];
    // int measure = 0;
    // for (int ii = 0; ii < 2 * Nsweep; ii++) {
    //   data[ii] = 0.0;
    // }

    // char filename[80];
    // sprintf(filename, "out_sweeptest.dat");
    // FILE *file = fopen(filename, "w");

    // int Nsw = Nsweep * N - 10;
    for (int i = 0; i < NT; i++) {

#ifdef _USE_OMP
#pragma omp parallel for
#endif
      for (int ii = 0; ii < N; ++ii) {
        for (j = 0; j < walker; ++j) {
          spin[j * N + ii] = spini[j * N + ii];
        }
      }
#ifdef _USE_OMP
#pragma omp parallel for
#endif
      for (int ii = 0; ii < walker; ii++) {
        En[ii] = Ei[ii];
        Mf[ii] = Mi[ii];
      }

      average(En, walker, &en);
      average(Mf, walker, &mag);

      t1 = cclock();

      for (int it = 0; it < N * Nsweep; ++it) {

        buff = it % N;
        transpose(spin, spinT, walker, N);

        cblas_dgemv(CblasRowMajor, CblasNoTrans, walker, N, 1.0, spin, N,
                    J + buff * N, 1, 0.0, dE, 1);
        buffer = spinT + buff * walker;

#ifdef _USE_OMP
#pragma omp parallel for
#endif
        for (int kk = 0; kk < walker; kk++) {
          dE[kk] = -4.0 * buffer[kk] * dE[kk];
        }

        metropolis(walker, N, buff, T, En, Mf, dE, spin);

        average(En, walker, &en);
        average(Mf, walker, &mag);
        // if (buff == N - 1) {
        //   printf("%f\t%f\t%f\n", T, en, mag);
        //   data[measure] = en;
        //   measure = measure + 1;
        //   data[measure] = mag;
        //   measure = measure + 1;
        // }
      }
      // measure = 0;
      // stder(T, data, 2 * Nsweep, file);
      T = T + dT;
    }
    t2 = cclock();
    printf("%d %f\n", N, (t2 - t1) / Nsweep);
    free(spin);
    free(spini);
    free(spinT);
    free(tmp);
  }
  // fclose(file);

  free(J);
  free(dE);
  free(En);
  free(Mi);
  free(Mf);
  return 0;
}
