#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #include "Matrixs.h"
#include "LU-factorisation.h"

#define N 10
#define MAXNUMBER 10
#define Mcount 13
#define DEBUG 0

typedef struct {
    double condition, error, rel_err, nev;
}conditionalDep;

void fullFree(Matrix** matrixs, int len) {
    for (int i = len-1; i >= 0; i--) {
        if (matrixs[i] != NULL) {
            freeMatrix(matrixs[i]);
        }
    }
}

conditionalDep ConditionalDepend(double cond, int dimension, FILE* output) {
    double n = dimension;
    long IterCount = 0;
    Matrix* matrixs[Mcount] = { 0 };

    Matrix* xfix = createXVector(n);
    matrixs[0] = xfix;
    printf("x*");
    printMatrix(xfix);

    Matrix *D = builtDMatrix(n, cond);
    matrixs[1] = D;
    printf("D");
    printMatrix(D);

    Matrix *Q = createQMatrix(n);
    matrixs[2] = Q;
    printf("Q");
    printMatrix(Q);

    Matrix* QD = multiplyMatrix(Q, D);
    matrixs[3] = QD;
    printf("QD");
    printMatrix(QD);

    Matrix* Qt = transposeMatrix(Q);
    matrixs[4] = Qt;
    printf("Qt");
    printMatrix(Qt);

    Matrix* A = multiplyMatrix(QD, Qt);
    matrixs[5] = A;
    printf("A");
    printMatrix(A);

    Matrix* b = multiplyMatrix(A, xfix);
    matrixs[6] = b;
    printf("b");
    printMatrix(b);

    Matrix* Comby = LU_factorisation(A, &IterCount); // LU = A
    matrixs[7] = Comby;
    printf("Comby");
    printMatrix(Comby);

    Matrix* y = forwardSM(Comby, b, &IterCount);
    matrixs[8] = y;
    printf("y");
    printMatrix(y);

    Matrix* x = backwardSM(Comby, y, &IterCount);
    matrixs[9] = x;
    printf("x");
    printMatrix(x);

    Matrix* errVector = addMatrix(x, xfix, '-');
    printf("err");
    printMatrix(errVector);
    matrixs[10] = errVector;
    double err = normMatrix(errVector);
    double err_norm = err / normMatrix(xfix);

    Matrix* b2 = multiplyMatrix(A, x);
    matrixs[11] = b2;
    Matrix* nevVector = addMatrix(b2, b, '-');
    matrixs[12] = nevVector;
    double nev = normMatrix(nevVector);
    double nev_norm = nev / normMatrix(b);

    fullFree(matrixs, Mcount);
    printf("cond\terr\t\t\t\terrd\t\t\tnev\tnev_norm\n");
    fprintf(output, "%e,%e,%e,%e,%e\n",cond , err, err_norm, nev, nev_norm);
    printf("%.0e\t% e\t% e\t% e\t% e\n",cond , err, err_norm, nev, nev_norm);

}

void DimensionDepend(double cond, int dimension, FILE* output) {
    long IterCount = 0;
    double n = dimension;
    Matrix* matrixs[Mcount] = { 0 };

    Matrix* xfix = createXVector(n);
    matrixs[0] = xfix;

    Matrix *D = builtDMatrix(n, cond);
    matrixs[1] = D;

    Matrix *Q = createQMatrix(n);
    matrixs[2] = Q;

    Matrix* QD = multiplyMatrix(Q, D);
    matrixs[3] = QD;

    Matrix* Qt = transposeMatrix(Q);
    matrixs[4] = Qt;

    Matrix* A = multiplyMatrix(QD, Qt);
    matrixs[5] = A;

    Matrix* b = multiplyMatrix(A, xfix);
    matrixs[6] = b;

    Matrix* Comby = LU_factorisation(A, &IterCount); // LU = A
    matrixs[7] = Comby;

    Matrix* y = forwardSM(Comby, b, &IterCount);
    matrixs[8] = y;

    Matrix* x = backwardSM(Comby, y, &IterCount);
    matrixs[9] = x;

    printf("cond\tdim\titer\n");
    printf("%.0e\t%d\t%ld\n",cond , dimension, IterCount);
    fprintf(output, "%e,%d,%ld\n", cond , dimension, IterCount);

    fullFree(matrixs, Mcount);

}


int main(void) {
    printf("Hello, World!\n");

    double cond = 1e1;
    int n = 10;

    // CONDITION DEPEND
    FILE* out = fopen("ConditionalDepend.csv", "w");
    if (out == NULL) {
        perror("ConditionalDepend.csv");
        return 1;
    }

    for (int p = 1; p < 16; p++) {
        double cond = pow(10,p);
        ConditionalDepend(cond, n, out);
    }
    fclose(out);

    // DIMENSION DEPEND
    FILE* out2 = fopen("DimensionDepend.csv", "w");
    if (out2 == NULL) {
        perror("DimensionDepend.csv");
        return 1;
    }

    for (int p = 0; p < 6; p++) {
        int dim = 10 * (2<<(p));
        // printf("%d\n", dim);
        DimensionDepend(cond, dim, out2);
    }
    fclose(out2);

    return 0;
}
