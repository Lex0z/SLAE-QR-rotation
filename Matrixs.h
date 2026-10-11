//
// Created by Алексей on 11.10.2026.
//

#ifndef LAB2_MATRIXS_H
#define LAB2_MATRIXS_H
#endif


typedef struct {
    int n; // row
    int m; // column
    double* d;
}Matrix;

Matrix* newMatrix(int n, int m);

Matrix* createMatrix(int n, int m, double* d);

void printMatrix(Matrix* matrix);

void freeMatrix(Matrix* matrix);

double normMatrix(Matrix* T);

Matrix* multiplyMatrix(Matrix* A, Matrix* B);

Matrix* createRandomMatrix(int n, int m);

Matrix* transposeMatrix(Matrix* M);

Matrix* addMatrix(Matrix* A, Matrix* B, char operation);
