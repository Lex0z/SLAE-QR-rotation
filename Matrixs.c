//
// Created by Алексей on 11.10.2026.
//

#include "Matrixs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DEBUG 0

Matrix* newMatrix(int n, int m) {
    Matrix *matrix = (Matrix*)malloc(sizeof(Matrix));
    matrix->n = n;
    matrix->m = m;
    matrix->d = (double *)malloc(sizeof(double) * m * n);
    return matrix;
}

Matrix* createMatrix(int n, int m, double* d) {
    Matrix *matrix = newMatrix(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (DEBUG) {printf("$%lf ", d[i * m + j]);}
            matrix->d[i * m + j] = d[i * m + j];
        }
    }
    return matrix;
}

void printMatrix(Matrix* matrix) {
    printf("\nMatrix %d x %d:\n", matrix->n, matrix->m);
    for (int i = 0; i < matrix->n; i++) {
        for (int j = 0; j < matrix->m; j++) {
            printf("% e ",matrix->d[i * matrix->m  + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeMatrix(Matrix* matrix) {
    free(matrix->d);
    free(matrix);
}

double normMatrix(Matrix* T) {
    double mx = 0.0;
    for (int j = 0; j < T->m; j++) {
        double sum = 0.0;
        for (int i = 0; i < T->n; i++) {
            sum += fabs(T->d[i*T->m + j]);
        }
        if (sum > mx) { mx = sum; }
    }
    return mx;
}

Matrix* multiplyMatrix(Matrix* A, Matrix* B) {
    int n = A->n;
    int p = A->m;
    int m = B->m;
    Matrix *C = newMatrix(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < p; k++) {
                sum += A->d[i * A->m + k] * B->d[k * B->m + j];
            }
            C->d[i * C->m + j] = sum;
        }
    }
    return C;
}

Matrix* createRandomMatrix(int n, int m) {
    double *d = (double *)malloc(sizeof(double) * n * m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            d[i*m + j] = (double) rand() / RAND_MAX;
        }
    }
    return createMatrix(n, m, d);
}

Matrix* transposeMatrix(Matrix* M) {
    double *d = (double *)malloc(sizeof(double) * M->n * M->m);
    for (int i = 0; i < M->n; i++) {
        for (int j = 0; j < M->m; j++) {
            d[j*M->n + i] = M->d[i*M->m + j];
        }
    }
    return createMatrix(M->m, M->n, d);
}

Matrix* addMatrix(Matrix* A, Matrix* B, char op) {
    if (A->m != B->m || A->n != B->n) {
        printf("Matrix dimensions do not match");
        return NULL;
    }
    double k = 0;
    switch (op) {
        case '+':
            k = 1;
            break;
        case '-':
            k = -1;
            break;
        default:
            printf("Invalid operation");
            return NULL;
            break;
    }

    Matrix* C = newMatrix(A->n, A->m);

    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->m; j++) {
            C->d[i*C->m + j] = A->d[i*A->m + j] + k * B->d[i*B->m + j];
        }
    }
    return C;
}