//
// Created by Алексей on 25.10.2026.
//

#include "QR.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



Matrix* backwardSM(Matrix* T, Matrix* y, long* IterCount) {
#define U(i, j) getR(T, i, j)
    int n = T->n;
    Matrix* x = newMatrix(n, 1);

    for (int i = n-1; i >= 0; --i) {
        double s = 0;
        for (int j = i+1; j < n; ++j) {
            s+= U(i,j)*x->d[j];
            (*IterCount)++;
        }
        x->d[i] = (y->d[i] - s)/U(i,i);
        (*IterCount)++;
    }

#undef U
    return x;
}

Matrix* createQMatrix(int n) {
    int m = 1;
    Matrix *W = createRandomMatrix(n,m);

    double norm2 = 0.0;
    for (int i = 0; i < n; i++) norm2 += W->d[i] * W->d[i];
    norm2 = sqrt(norm2);
    for (int i = 0; i < n; i++) W->d[i] /= norm2;

    m = n;
    Matrix *Q = newMatrix(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double e = 0;
            if (i == j) {e = 1;}
            Q->d[i * m + j] = e - 2 * W->d[i] * W->d[j];
        }
    }

    freeMatrix(W);
    return Q;
}

Matrix* builtDMatrix(int n, double cond) {
    double* d = (double *)malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        d[i] = 1 + (cond - 1) * i / (n - 1);
    }
    int m = n;
    Matrix* D = newMatrix(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (i == j) {D->d[i * m + j] = d[i];}
            else { D->d[i * m + j] = 0; }
        }
    }

    free(d);
    return D;
}

Matrix* createXVector(int n) {
    double* d = (double *)malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        d[i] = (double)i+1;
    }

    Matrix* x = createMatrix(n, 1, d);
    free(d);
    return x;
}