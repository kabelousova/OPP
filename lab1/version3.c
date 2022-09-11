#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>

#define eps 1E-5
#define tau 1E-4

#define NORM 0
#define VECTOR 1

#define DOUBLE_RANGE 8

double *genA(int N) {
    double *A = malloc(sizeof(double) * N * N);

    srand(time(NULL));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            A[i * N + j] = (((double) rand()) / RAND_MAX) * 2 * DOUBLE_RANGE - DOUBLE_RANGE;
            A[j * N + i] = A[i * N + j];
        }

        A[i * N + i] += 300;
    }

    return A;
}

double *genB(int N) {
    double *B = malloc(sizeof(double) * N);

    for (int i = 0; i < N; ++i) {
        B[i] = 0.0001;
    }

    return B;
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

void distributionVector(double *v, double *rv, const int *countColumns) {
    int size, rank;
    int *countData, *shifts;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    countData = (int *) malloc(sizeof(int) * size);
    shifts = (int *) malloc(sizeof(int) * size);

    for (int i = 0; i < size; ++i) {
        countData[rank] = countColumns[rank];
    }

    shifts[0] = 0;
    for (int i = 1; i < size; ++i) {
        shifts[rank] = shifts[rank - 1] + countData[rank - 1];
    }

    MPI_Scatterv(v, countData, shifts, MPI_DOUBLE, rv, countData[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(shifts);
    free(countData);
}

// Главные процесс (0) отправляет остальным их куски матрицы и вектора b
void distribution(int N, double *A, double *rA, double *b, double *rb, const int *countColumns) {
    int size, rank;
    int *countData, *shifts;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    countData = (int *) malloc(sizeof(int) * size);
    shifts = (int *) malloc(sizeof(int) * size);

    for (int i = 0; i < size; ++i) {
        countData[rank] = countColumns[rank];
    }

    shifts[0] = 0;
    for (int i = 1; i < size; ++i) {
        shifts[rank] = shifts[rank - 1] + countData[rank - 1];
    }

    MPI_Scatterv(b, countData, shifts, MPI_DOUBLE, rb, countData[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD); //для вектора b

    for (int i = 0; i < size; ++i) {
        countData[rank] *= N;
        shifts[rank] *= N;
    }

    MPI_Scatterv(A, countData, shifts, MPI_DOUBLE, rA, countData[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD); //для матрицы

    free(shifts);
    free(countData);
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

    double *A = NULL, *rA = NULL;
    double *x = NULL;
    double *b = NULL, *rb = NULL;
    double *c = malloc(sizeof(double) * N);
    double *completeC = malloc(sizeof(double) * N);

    // Массив кол-ва столбцов по процессам
    int *countColumns = malloc(sizeof(int) * size);

    double smallNorm, normb = 0;

    // Номер первого столбца в нашем процессе
    int shift = 0;
    int remain = N % size;

    // Определяем, сколько столбцов у каждого процесса
    for (int i = 0; i < size; ++i) {
        countColumns[i] = i < remain ? N / size + 1 : N / size;
        if (i < rank) {
            shift += countColumns[i];
        }
    }

    if (rank == 0) {
        A = genA(N);
        b = genB(N);
    }

    rA = malloc(sizeof(double) * N * countColumns[rank]);
    rb = malloc(sizeof(double) * countColumns[rank]);

    distribution(N, A, rA, b, rb, countColumns);

    x = malloc(sizeof(double) * countColumns[rank]);

    for (int i = 0; i < countColumns[rank]; ++i) {
        x[i] = 0;
    }

    // Cчитаем квадраты кусков b
    smallNorm = sq_sum(countColumns[rank], rb);

    // Собираем квадраты в одну переменную
    MPI_Allreduce(&smallNorm, &normb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Так как сейчас это квадрат норм, то извлекаем корень
    normb = sqrt(normb);

    while (count < 5) {
        double bignorm;
        mul(N, countColumns[rank], rA, x, c);

        // Собрали результат перемножения
        memset(completeC, 0, N * sizeof(
        double *N));
        MPI_Allreduce(c, completeC, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // Разослали всем их куски
        distributionVector(completeC, c, countColumns);

        sub(countColumns[rank], c, rb);

        smallNorm = sq_sum(countColumns[rank], c);
        // Собираем квадраты в одну переменную
        MPI_Allreduce(&smallNorm, &bignorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        bignorm = sqrt(bignorm);

        mulConst(countColumns[rank], c, tau);

        sub(countColumns[rank], x, c);

        // Считаем условие выхода
        if (bignorm / normb < eps) {
            count += 1;
        } else count = 0;
    }

    MPI_Finalize();

    free(countColumns);
    free(A);
    free(rA);
    free(x);
    free(b);
    free(rb);
    free(c);
    free(completeC);

    return 0;
}
