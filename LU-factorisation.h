//
// Created by Алексей on 11.10.2026.
//

#ifndef LAB2_LU_FACTORISATION_H
#define LAB2_LU_FACTORISATION_H

#endif //LAB2_LU_FACTORISATION_H

#include "Matrixs.h"

Matrix* LU_factorisation(Matrix *A, long* IterCount);

double getL(Matrix* T, int i , int j);

double getR(Matrix* T, int i , int j);

Matrix* forwardSM(Matrix* T, Matrix* b, long* IterCount); //forward substitution method

Matrix* backwardSM(Matrix* T, Matrix* b, long* IterCount); //backward substitution method

Matrix* createQMatrix(int n);

Matrix* builtDMatrix(int n, double cond);

Matrix* createXVector(int n);