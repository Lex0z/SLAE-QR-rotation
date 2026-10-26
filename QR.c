#include "QR.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

static Matrix* copyMatrix(const Matrix* source) {
    if (source == NULL || source->n <= 0 || source->m <= 0 || source->d == NULL) {
        return NULL;
    }
    Matrix* copy = newMatrix(source->n, source->m);
    if (copy == NULL || copy->d == NULL) {
        if (copy != NULL) freeMatrix(copy);
        return NULL;
    }
    for (int i = 0; i < source->n * source->m; ++i) {
        copy->d[i] = source->d[i];
    }
    return copy;
}

Matrix* backwardSM(const Matrix* T, const Matrix* y, long* IterCount) {
    if (T == NULL || y == NULL || T->n != T->m || y->n != T->n || y->m != 1) {
        return NULL;
    }

    const int n = T->n;
    Matrix* x = newMatrix(n, 1);
    if (x == NULL || x->d == NULL) {
        if (x != NULL) freeMatrix(x);
        return NULL;
    }

    double scale = 0.0;
    for (int i = 0; i < n; ++i) {
        double value = fabs(T->d[i * T->m + i]);
        if (value > scale) scale = value;
    }
    const double tolerance = DBL_EPSILON * (double)n * scale;

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += T->d[i * T->m + j] * x->d[j];
            if (IterCount != NULL) *IterCount += 2; // multiply and addition
        }

        const double diagonal = T->d[i * T->m + i];
        if (!isfinite(diagonal) || fabs(diagonal) <= tolerance) {
            freeMatrix(x);
            return NULL;
        }
        x->d[i] = (y->d[i] - sum) / diagonal;
        if (IterCount != NULL) *IterCount += 2; // subtraction and division
    }
    return x;
}

Matrix* QRsolve(const Matrix* A, const Matrix* b, long* IterCount) {
    if (A == NULL || b == NULL || A->n <= 0 || A->n != A->m ||
        b->n != A->n || b->m != 1 || A->d == NULL || b->d == NULL) {
        return NULL;
    }

    Matrix* T = copyMatrix(A);
    Matrix* y = copyMatrix(b);
    if (T == NULL || y == NULL) {
        if (T != NULL) freeMatrix(T);
        if (y != NULL) freeMatrix(y);
        return NULL;
    }

    const int n = T->n;
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            const double a = T->d[k * n + k];
            const double d = T->d[i * n + k];
            if (d == 0.0) continue;

            const double rho = hypot(a, d);
            if (rho == 0.0) continue;
            const double c = a / rho;
            const double s = d / rho;
            if (IterCount != NULL) *IterCount += 3; // hypot, two divisions

            for (int j = k; j < n; ++j) {
                const double upper = T->d[k * n + j];
                const double lower = T->d[i * n + j];
                T->d[k * n + j] = c * upper + s * lower;
                T->d[i * n + j] = -s * upper + c * lower;
                if (IterCount != NULL) *IterCount += 6; // four products, two additions
            }

            const double upperB = y->d[k];
            const double lowerB = y->d[i];
            y->d[k] = c * upperB + s * lowerB;
            y->d[i] = -s * upperB + c * lowerB;
            if (IterCount != NULL) *IterCount += 6;

            /* Avoid retaining a round-off remnant in the eliminated position. */
            T->d[i * n + k] = 0.0;
        }
    }

    Matrix* x = backwardSM(T, y, IterCount);
    freeMatrix(T);
    freeMatrix(y);
    return x;
}

Matrix* createQMatrix(int n) {
    if (n <= 0) return NULL;
    Matrix* W = createRandomMatrix(n, 1);
    if (W == NULL) return NULL;

    double norm2 = 0.0;
    for (int i = 0; i < n; ++i) norm2 += W->d[i] * W->d[i];
    norm2 = sqrt(norm2);
    if (norm2 == 0.0) {
        freeMatrix(W);
        return NULL;
    }
    for (int i = 0; i < n; ++i) W->d[i] /= norm2;

    Matrix* Q = newMatrix(n, n);
    if (Q == NULL || Q->d == NULL) {
        if (Q != NULL) freeMatrix(Q);
        freeMatrix(W);
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double identity = (i == j) ? 1.0 : 0.0;
            Q->d[i * n + j] = identity - 2.0 * W->d[i] * W->d[j];
        }
    }
    freeMatrix(W);
    return Q;
}

Matrix* builtDMatrix(int n, double cond) {
    if (n <= 0 || !isfinite(cond) || cond < 1.0) return NULL;
    Matrix* D = newMatrix(n, n);
    if (D == NULL || D->d == NULL) {
        if (D != NULL) freeMatrix(D);
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) D->d[i * n + j] = 0.0;
        D->d[i * n + i] = 1.0 + (cond - 1.0) * (double)i / (double)(n - 1);
    }
    return D;
}

Matrix* createXVector(int n) {
    if (n <= 0) return NULL;
    Matrix* x = newMatrix(n, 1);
    if (x == NULL || x->d == NULL) {
        if (x != NULL) freeMatrix(x);
        return NULL;
    }
    for (int i = 0; i < n; ++i) x->d[i] = (double)i + 1.0;
    return x;
}
