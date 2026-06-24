#ifndef ARCTAN2_APPROX_H
#define ARCTAN2_APPROX_H

#include "matrix_utils.h"

static constexpr int CORDIC_ITERATIONS = 24;

double manual_sqrt(double x);
double cordic_atan(double y, double x);
double manual_atan2(double y, double x);
double manual_pi();

Matrix h_joint(const Matrix& x_joint);
Matrix H_jacobian_joint(const Matrix& x_joint);
Matrix h_full(const Matrix& x_full);
Matrix H_jacobian_full(const Matrix& x_full);

#endif // ARCTAN2_APPROX_H
