#include <iostream>
#include <mpi.h>
#include <malloc.h>

#define COLUMNS 2
#define ROWS 8
#define N 2080
#define M 2080
#define K 2080
#define countRows (N / ROWS)
#define countColumns (K / COLUMNS)

using namespace
std;

void sendA(double **A, int rank, MPI_Comm line, MPI_Comm col, MPI_Comm comm) {
    int *coord = new
    int[2]; //одномерный массив для координат
    double *tA = (double *) malloc(sizeof(double) * M * countRows);

    MPI_Cart_coords(comm, rank, 2, coord); //получаем координаты текущего

    if (!coord[1]) {
        MPI_Scatter(*A, M * countRows, MPI_DOUBLE, tA, M * countRows, MPI_DOUBLE, 0, col); //рассылаем как на рисунке 1
    }

    free(*A);
    *A = tA;

    MPI_Bcast(*A, M * countRows, MPI_DOUBLE, 0, line); //рассылка по строкам
}

void sendB(double **B, int rank, MPI_Comm line, MPI_Comm colu, MPI_Comm comm) {
    int *coord = new
    int[2]; //одномерный массив для координат
    double *tB = (double *) malloc(sizeof(double) * M * countColumns);
    int countb[COLUMNS]; // сколько элементов типа col отправить процессам
    int dispb[COLUMNS]; //смещение каждого блока относительно начала

    MPI_Cart_coords(comm, rank, 2, coord);

    MPI_Datatype col, types[2];
    MPI_Type_vector(M, countColumns, K, MPI_DOUBLE, &types[0]);

    MPI_Aint sizeofdouble, disp[2], lb;
    int blen[2];
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &sizeofdouble); //получаем информацию о типе данных
    blen[0] = 1;
    blen[1] = 1;
    disp[0] = 0;
    disp[1] = sizeofdouble * countColumns;
    types[1] = MPI_UB;
    MPI_Type_create_struct(2, blen, disp, types, &col); //создаем структуру
    MPI_Type_commit(&col); //регистрируем созданный тип данных col

    for (int i = 0; i < COLUMNS; i++) {
        dispb[i] = i;
        countb[i] = 1;
    }

    if (!coord[0]) {
        MPI_Scatterv(*B, countb, dispb, col, tB, countColumns * M, MPI_DOUBLE, 0, line);
    } //раздается как на рисунке 2

    free(*B);
    *B = tB;

    MPI_Bcast(*B, M * countColumns, MPI_DOUBLE, 0, colu); //рассылка уже по столбцам

    MPI_Type_free(&col);
    MPI_Type_free(types); //удаляем типы
}

void mulC(double *A, double *B, double **C) {

    *C = (double *) malloc(sizeof(double) * countRows * countColumns);

    for (int i = 0; i < countRows; ++i) {
        for (int j = 0; j < countColumns; ++j) {
            (*C)[i * countColumns + j] = 0.0;
            for (int k = 0; k < M; ++k) {
                (*C)[i * countColumns + j] += A[i * M + k] * B[k * countColumns + j];
            }
        }
    }
}

void gatherC(double **C, MPI_Comm comm, int rank) {
    double *tC = !rank ? (double *) malloc(sizeof(double) * N * K) : NULL;
    int *dispc = (int *) malloc(sizeof(int) * ROWS * COLUMNS); //сдвиг данных от начала массива
    int *countc = (int *) malloc(sizeof(int) * ROWS * COLUMNS); //кол-во данных

    MPI_Datatype col, types[2];
    MPI_Type_vector(M, countColumns, K, MPI_DOUBLE, &types[0]);

    MPI_Aint sizeofdouble, disp[2], lb;
    int blen[2];
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &sizeofdouble);//double
    blen[0] = 1;
    blen[1] = 1;
    disp[0] = 0;
    disp[1] = sizeofdouble * countColumns;
    types[1] = MPI_UB;
    MPI_Type_create_struct(2, blen, disp, types, &col);
    MPI_Type_commit(&col); //регистрируем новый тип

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            dispc[i * COLUMNS + j] = (i * COLUMNS * countRows + j); //ном
            countc[i * COLUMNS + j] = 1;
        }
    }

    MPI_Gatherv(*C, countRows * countColumns, MPI_DOUBLE, tC, countc, dispc, col, 0, comm);
    free(*C);
    *C = tC;
}

void createAB(double **A, double **B, double **C, int rank) {
    if (!rank) {
        *A = (double *) malloc(sizeof(double) * N * M);
        *B = (double *) malloc(sizeof(double) * M * K);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                (*A)[M * i + j] = M * i + j;
            }
        }

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < K; ++j) {
                (*B)[K * i + j] = (i == j);
            }
        }
    }
}


int main(int argc, char **argv) {
    double *A = NULL;
    double *B = NULL;
    double *C = NULL;

    double start, end;

    int rank;
    int per[2] = {0, 0};
    int l[2] = {0, 1}; //строка
    int c[2] = {1, 0}; //столбец
    int m[2] = {ROWS, COLUMNS}; //размер решетки

    MPI_Init(&argc, &argv);

    MPI_Comm decart, column, row;
    MPI_Cart_create(MPI_COMM_WORLD, 2, m, per, 0, &decart);
    MPI_Cart_sub(decart, l, &row); //из решетки получить по строкам
    MPI_Cart_sub(decart, c, &column); //по столбцам

    MPI_Comm_rank(decart, &rank);

    createAB(&A, &B, &C, rank);
    start = MPI_Wtime();
    sendA(&A, rank, row, column, decart);
    sendB(&B, rank, row, column, decart);
    mulC(A, B, &C);
    gatherC(&C, decart, rank);

    end = MPI_Wtime();

    if (!rank) {
        cout << end - start << endl;
    }

    MPI_Finalize();
    return 0;
}
