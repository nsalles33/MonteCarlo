#ifndef MOD_H
#define MOD_H

#include <cblas.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
// #include <unistd.h>

double cclock();
void transpose(double *start, double *end, int N, int M);
void metropolis(int M, int N, int step, double T, double *en, double *mag,
                double *dE, double *sigma);
void average(double *x,int N, double *out);
void magnetization(double *sigma, double *mag, int N, int idx);
void stder(double T,double *data, int N, FILE *file);
#endif
