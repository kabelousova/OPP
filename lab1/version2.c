#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#define eps 1E-5
#define tau 1E-4

#define NORM 0
#define VECTOR 1

#define DOUBLE_RANGE 8

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

void mul(int M, int N, const double *A, const double *x, double *c) {
    for (int i = 0; i < M; i++) {
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

double sq_sum(int N, const double *v) {
    double norm = 0;

    for (int i = 0; i < N; i++) {
        norm += v[i] * v[i];
    }

    return norm;
}

void rootProcess(int *count, double *x, const int *countRows, double bignorm, double normb) {
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 1, shift = countRows[rank]; i < size; shift += countRows[i], ++i) {
        double tmp;

        MPI_Recv(&x[shift], countRows[i], MPI_DOUBLE, i, VECTOR, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Recv(&tmp, 1, MPI_DOUBLE, i, NORM, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        bignorm += tmp;
    }

    bignorm = sqrt(bignorm);

    printf("%lf %lf %lf\n", bignorm / normb, bignorm, normb);

    if (bignorm / normb < eps) {
        *count += 1;
    } else *count = 0;
}

// Главные процесс (0) отправляет остальным их куски матрицы
void distribution(int N, double **A, const int *countRows) {
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank) {
        *A = genA(N);

        for (int i = 1, shift = countRows[rank]; i < size; shift += countRows[i], ++i) {
            MPI_Send(*A + shift * N, countRows[i] * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    } else {
        *A = malloc(sizeof(double) * countRows[rank] * N);
        MPI_Recv(*A, countRows[rank] * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    }
}

int main(int argc, char **argv) {
    int N = 1000;
    int count = 0;

    // Номер нашего процесса и их кол-во
    int rank, size;

    // Инициализираум MPI и получаем номер процесса и их кол-во
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *A = NULL;
    double *x = malloc(sizeof(double) * N);
    double *b = malloc(sizeof(double) * N);
    double *c = NULL;

    // Массив кол-ва строк по процессам
    int *countRows = malloc(sizeof(int) * size);

    double normb;

    // Номер первой строки в нашем процессе
    int shift = 0;
    int remain = N % size;

    // Определяем, сколько строк у каждого процесса
    for (int i = 0; i < size; ++i) {
        countRows[i] = i < remain ? N / size + 1 : N / size;
        if (i < rank) {
            shift += countRows[i];
        }
    }

    distribution(N, &A, countRows);

    c = malloc(sizeof(double) * countRows[rank]);

    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        b[i] = 0.0001;
    }

    normb = sqrt(sq_sum(N, b));

    while (count < 5) {
        double bignorm;
        mul(countRows[rank], N, A, x, c);
        sub(countRows[rank], c, b + shift);

        bignorm = sq_sum(countRows[rank], c);

        mulConst(countRows[rank], c, tau);

        sub(countRows[rank], x + shift, c);

        // Теперь отправляем посчитанный кусок вектора и норму куска нулевому процессу
        if (rank) {
            MPI_Send(&bignorm, 1, MPI_DOUBLE, 0, NORM, MPI_COMM_WORLD);
            MPI_Send(x + shift, countRows[rank], MPI_DOUBLE, 0, VECTOR, MPI_COMM_WORLD);
        } else {
            rootProcess(&count, x, countRows, bignorm, normb);
        }

        // Групповая рассылка всем count и собранного в нулевом процессе x
        MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    free(countRows);
    free(A);
    free(x);
    free(b);
    free(c);

    return 0;
}
