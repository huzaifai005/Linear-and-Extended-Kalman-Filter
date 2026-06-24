#ifndef KALMAN_MATRICES_H
#define KALMAN_MATRICES_H

#include "matrix_utils.h"

static constexpr int N_AXIS = 4;
static constexpr int N_JOINT = 12;
static constexpr int N_JOINTS = 23;
static constexpr int N_STATE = 276;
static constexpr int N_MEAS = 69;

Matrix make_F1(double dt);
Matrix make_Q1(double dt, double phi);

Matrix make_F_joint(double dt);
Matrix make_Q_joint(double dt, double phi);
Matrix make_H_joint();

Matrix make_F_full(double dt);
Matrix make_Q_full(double dt, double phi);
Matrix make_H_full();
Matrix make_R_full(double r_noise);
Matrix make_P0_full(double p_pos = 1.0, double p_vel = 100.0,
double p_acc = 100.0, double p_jerk = 100.0);
Matrix make_x0(const std::vector<double>& z0);

#endif
