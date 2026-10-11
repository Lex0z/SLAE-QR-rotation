//
// Created by Алексей on 11.10.2026.
//

#include "LU-factorisation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

Matrix* LU_factorisation(Matrix* A, long* IterCount) {
    if (A->n != A->m) {
        printf("Error in LU_factorisation: n != m");
        return NULL;
    }
    int n = A->n;
    Matrix* T = newMatrix(n,n);
    for (int m = 0; m < n; m++) {
        for (int j = m; j < n; j++) {
            double s = 0;
            for (int k = 0; k < m; ++k) {
                s+=getL(T,m, k) * getR(T, k, j);
                (*IterCount)++;
            }
            T->d[m * T->m + j] = A->d[m * A->m + j] - s; // U_{mj}
            (*IterCount)++;
            if (m == j && T->d[m * T->m + j] == 0.0 ) {
                printf("Error in LU_factorisation: U_{%d%d} = 0", m,j);
                return NULL;
            }
        }

        for (int i = m+1; i < n; ++i) {
            double s = 0;
            for (int k = 0; k < m; ++k) {
                s+=getL(T,i,k) * getR(T,k,m);
                (*IterCount)++;
            }
            T->d[i*n + m] = (A->d[i*n + m] - s) / getR(T,m,m); // L_{im}
            (*IterCount)++;
        }
    }
    return T;
}

double getL(Matrix* T, int i , int j) {
    if (i == j) {
        return 1.0;
    }
    if (i>j) {
        return T->d[i*T->m + j];
    } else {
        return 0.0;
    }
}

double getR(Matrix* T, int i , int j) {
    if (i>j) {
        return 0.0;
    } else {
        return T->d[i*T->m + j];
    }
}

Matrix* forwardSM(Matrix* T, Matrix* b, long* IterCount) {
    #define L(i, j) getL(T, i, j)

    int n = T->n;
    Matrix* y = newMatrix(n,1);

    for (int i = 0; i < n; ++i) {
        double s = 0;
        for (int j = 0; j < i; ++j) {
            s+= L(i,j)*y->d[j];
            (*IterCount)++;
        }
        y->d[i] = b->d[i] - s;
        (*IterCount)++;
    }

    #undef L
    return y;
}

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

Matrix* builtDmatrix(int n, double cond) {
    double* d = (double *)malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        d[i] = 1 + (cond - 1) * i / n;
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
    return createMatrix(n, 1, d);
}
