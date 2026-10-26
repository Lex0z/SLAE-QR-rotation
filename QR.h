//
// Created by Алексей on 25.10.2026.
//

#ifndef SLAE_LU_METHOD_QH_H
#define SLAE_LU_METHOD_QH_H

#endif //SLAE_LU_METHOD_QH_H

#include "Matrixs.h"



Matrix* backwardSM(Matrix* T, Matrix* b, long* IterCount); //backward substitution method

Matrix* createQMatrix(int n);

Matrix* builtDMatrix(int n, double cond);

Matrix* createXVector(int n);