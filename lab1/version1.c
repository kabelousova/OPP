#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define eps 1E-5
#define tau 1E-4

#define NORM 0
#define VECTOR 1

#define DOUBLE_RANGE 8 //отвечает за генерацию матрицы, в каком диапазоне будет

double *genA(int N) {
    double *A = malloc(sizeof(double) * N * N);

    srand(time(NULL));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            A[i * N + j] = (((double) rand()) / RAND_MAX) * 2 * DOUBLE_RANGE - DOUBLE_RANGE;
            A[j * N + i] = A[i * N + j];
        }

        A[i * N + i] += 300;
    }

    return A;
}

void mul(int N, const double *A, const double *x, double *c) {
    for (int i = 0; i < N; i++) {
        c[i] = 0;

        for (int j = 0; j < N; j++) {
            c[i] += A[i * N + j] * x[j];
        }
    }
}

void sub(int N, double *c, const double *b) {
    for (int i = 0; i < N; i++) {
        c[i] -= b[i];
    }
}

void mulConst(int N, double *c, double t) {
    for (int i = 0; i < N; i++) {
        c[i] *= t;
    }
}

double norm(int N, const double *v) {
    double norm = 0;

    for (int i = 0; i < N; i++) {
        norm += v[i] * v[i];
    }

    return sqrt(norm);
}

int main(int argc, char **argv) {
    int N = 1000;
    int count = 0;

    double *A = genA(N);
    double *x = malloc(sizeof(double) * N);
    double *b = malloc(sizeof(double) * N);
    double *c = malloc(sizeof(double) * N);

    double normb;

    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        b[i] = 0.0001;
    }

    normb = norm(N, b);

    while (count < 5) {
        double bignorm;
        mul(N, N, A, x, c);
        sub(N, c, b);

        bignorm = norm(N, c);

        mulConst(N, c, tau);

        sub(N, x, c);

        if (bignorm / normb < eps) {
            count += 1;
        } else count = 0;
    }

    free(A);
    free(x);
    free(b);
    free(c);

    return 0;
}
