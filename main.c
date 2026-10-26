#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "QR.h"

#ifndef PROJECT_ROOT_DIR
#define PROJECT_ROOT_DIR "."
#endif

static Matrix* buildAMatrix(int n, double cond) {
    Matrix* Q = createQMatrix(n);
    Matrix* D = builtDMatrix(n, cond);
    if (Q == NULL || D == NULL) {
        freeMatrix(Q);
        freeMatrix(D);
        return NULL;
    }

    Matrix* QD = multiplyMatrix(Q, D);
    Matrix* Qt = transposeMatrix(Q);
    Matrix* A = NULL;
    if (QD != NULL && Qt != NULL) A = multiplyMatrix(QD, Qt);

    freeMatrix(Q);
    freeMatrix(D);
    freeMatrix(QD);
    freeMatrix(Qt);
    return A;
}

static int ConditionalDepend(double cond, int n, FILE* output) {
    long operationCount = 0;
    /* Keep Q identical across the condition-number sweep. */
    srand(1);
    Matrix* A = buildAMatrix(n, cond);
    Matrix* xExact = createXVector(n);
    Matrix* b = (A != NULL && xExact != NULL) ? multiplyMatrix(A, xExact) : NULL;
    Matrix* x = (A != NULL && b != NULL) ? QRsolve(A, b, &operationCount) : NULL;
    Matrix* error = (x != NULL && xExact != NULL) ? addMatrix(x, xExact, '-') : NULL;
    Matrix* Ax = (A != NULL && x != NULL) ? multiplyMatrix(A, x) : NULL;
    Matrix* residual = (Ax != NULL && b != NULL) ? addMatrix(Ax, b, '-') : NULL;

    if (A == NULL || xExact == NULL || b == NULL || x == NULL || error == NULL || Ax == NULL || residual == NULL) {
        fprintf(stderr, "Ошибка при решении системы: n=%d, cond=%g\n", n, cond);
        freeMatrix(A); freeMatrix(xExact); freeMatrix(b); freeMatrix(x);
        freeMatrix(error); freeMatrix(Ax); freeMatrix(residual);
        return 0;
    }

    const double errorNorm = normVector2(error);
    const double exactNorm = normVector2(xExact);
    const double residualNorm = normVector2(residual);
    const double bNorm = normVector2(b);
    const double relativeError = (exactNorm > 0.0) ? errorNorm / exactNorm : NAN;
    const double relativeResidual = (bNorm > 0.0) ? residualNorm / bNorm : residualNorm;

    /* CSV columns: cond, absolute error, relative error, residual, relative residual */
    fprintf(output, "%.12e,%.12e,%.12e,%.12e,%.12e\n",
            cond, errorNorm, relativeError, residualNorm, relativeResidual);
    printf("cond=%g, n=%d, rel_error=%.4e, rel_residual=%.4e\n",
           cond, n, relativeError, relativeResidual);

    freeMatrix(A); freeMatrix(xExact); freeMatrix(b); freeMatrix(x);
    freeMatrix(error); freeMatrix(Ax); freeMatrix(residual);
    return 1;
}

static int DimensionDepend(double cond, int n, FILE* output) {
    long operationCount = 0;
    Matrix* A = buildAMatrix(n, cond);
    Matrix* xExact = createXVector(n);
    Matrix* b = (A != NULL && xExact != NULL) ? multiplyMatrix(A, xExact) : NULL;
    Matrix* x = (A != NULL && b != NULL) ? QRsolve(A, b, &operationCount) : NULL;

    if (A == NULL || xExact == NULL || b == NULL || x == NULL) {
        fprintf(stderr, "Ошибка при решении системы: n=%d, cond=%g\n", n, cond);
        freeMatrix(A); freeMatrix(xExact); freeMatrix(b); freeMatrix(x);
        return 0;
    }

    /* CSV columns: cond, dimension, arithmetic operation count */
    fprintf(output, "%.12e,%d,%ld\n", cond, n, operationCount);
    printf("cond=%g, n=%d, operations=%ld\n", cond, n, operationCount);

    freeMatrix(A); freeMatrix(xExact); freeMatrix(b); freeMatrix(x);
    return 1;
}

static int openOutput(const char* relativePath, char* path, size_t pathSize, FILE** output) {
    snprintf(path, pathSize, "%s/%s", PROJECT_ROOT_DIR, relativePath);
    *output = fopen(path, "w");
    if (*output == NULL) {
        perror(path);
        return 0;
    }
    return 1;
}

int main(void) {
    char path[1024];
    FILE* output = NULL;
    const int conditionDimension = 10;
    const double dimensionCondition = 10.0;
    const size_t conditionCount = 15;

    srand(1); /* reproducible random Householder vectors */

    if (!openOutput("data/ConditionalDepend.csv", path, sizeof(path), &output)) return 1;
    for (size_t i = 0; i < conditionCount; ++i) {
        double cond = pow(10.0,(int)(i+1));
        if (!ConditionalDepend(cond, conditionDimension, output)) {
            fclose(output);
            return 1;
        }
    }
    fclose(output);

    if (!openOutput("data/DimensionDepend.csv", path, sizeof(path), &output)) return 1;
    for (int k = 1; k <= 20; ++k) {
        const int n = 15 * k;
        if (!DimensionDepend(dimensionCondition, n, output)) {
            fclose(output);
            return 1;
        }
    }
    fclose(output);
    return 0;
}
