#ifndef SLAE_QR_ROTATION_QR_H
#define SLAE_QR_ROTATION_QR_H

#include "Matrixs.h"

/* Solves A*x=b by Givens QR. A and b are not modified. Returns NULL
   if dimensions are invalid or the triangular system is numerically singular. */
Matrix* QRsolve(const Matrix* A, const Matrix* b, long* IterCount);

/* Back substitution for an upper-triangular matrix T and transformed RHS y. */
Matrix* backwardSM(const Matrix* T, const Matrix* y, long* IterCount);

/* Helpers for reproducible numerical experiments. */
Matrix* createQMatrix(int n);
Matrix* builtDMatrix(int n, double cond);
Matrix* createXVector(int n);

#endif
