// kalman_matrices.cpp
// Team Thunder — Kalman Filter Milestone 2

#include "kalman_matrices.h"
#include <cmath>
#include <stdexcept>

//   F1 = [ 1  dt  dt²/2  dt³/6 ]
//        [ 0   1    dt   dt²/2 ]
//        [ 0   0     1     dt  ]
//        [ 0   0     0      1  ]
Matrix make_F1(double dt) {
    Matrix F1(N_AXIS, N_AXIS);

    double dt2 = dt * dt;
    double dt3 = dt2 * dt;

    F1(0, 0) = 1.0;
    F1(0, 1) = dt;
    F1(0, 2) = dt2 / 2.0;
    F1(0, 3) = dt3 / 6.0;

    F1(1, 1) = 1.0;
    F1(1, 2) = dt;
    F1(1, 3) = dt2 / 2.0;

    F1(2, 2) = 1.0;
    F1(2, 3) = dt;

    F1(3, 3) = 1.0;

    return F1;
}

//   Q1 = Φ * [ dt⁷/252  dt⁶/72  dt⁵/30  dt⁴/24 ]
//            [ dt⁶/72   dt⁵/20  dt⁴/8   dt³/6  ]
//            [ dt⁵/30   dt⁴/8   dt³/3   dt²/2  ]
//            [ dt⁴/24   dt³/6   dt²/2   dt     ]
Matrix make_Q1(double dt, double phi) {
    
    Matrix Q1(N_AXIS, N_AXIS);

    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double dt5 = dt4 * dt;
    double dt6 = dt5 * dt;
    double dt7 = dt6 * dt;

    Q1(0, 0) = phi * dt7 / 252.0;
    Q1(0, 1) = phi * dt6 / 72.0;
    Q1(0, 2) = phi * dt5 / 30.0;
    Q1(0, 3) = phi * dt4 / 24.0;

    Q1(1, 0) = Q1(0, 1);
    Q1(1, 1) = phi * dt5 / 20.0;
    Q1(1, 2) = phi * dt4 / 8.0;
    Q1(1, 3) = phi * dt3 / 6.0;

    Q1(2, 0) = Q1(0, 2);
    Q1(2, 1) = Q1(1, 2);
    Q1(2, 2) = phi * dt3 / 3.0;
    Q1(2, 3) = phi * dt2 / 2.0;

    Q1(3, 0) = Q1(0, 3);
    Q1(3, 1) = Q1(1, 3);
    Q1(3, 2) = Q1(2, 3);
    Q1(3, 3) = phi * dt;

    return Q1;
}

// Per-joint matrices (12×12)
Matrix make_F_joint(double dt) {

    Matrix F1 = make_F1(dt);
    return block_diagonal({ F1, F1, F1 });
}

Matrix make_Q_joint(double dt, double phi) {

    Matrix Q1 = make_Q1(dt, phi);
    return block_diagonal({ Q1, Q1, Q1 });
}

Matrix make_H_joint() {
    
    // State layout: [px vx ax jx | py vy ay jy | pz vz az jz]
    //                0  1  2  3    4  5  6  7    8  9 10 11
   
    Matrix H(3, N_JOINT);
    H(0, 0) = 1.0;  
    H(1, 4) = 1.0;  
    H(2, 8) = 1.0;  
    return H;
}

// Full-body matrices (276×276 and 69×276)
Matrix make_F_full(double dt) {
    
    // Block diagonal: F_joint repeated 23 times
    Matrix F_joint = make_F_joint(dt);
    std::vector<Matrix> blocks(N_JOINTS, F_joint);
    return block_diagonal(blocks);
}

Matrix make_Q_full(double dt, double phi) {
    
    // Block diagonal: Q_joint repeated 23 times
    Matrix Q_joint = make_Q_joint(dt, phi);
    std::vector<Matrix> blocks(N_JOINTS, Q_joint);
    return block_diagonal(blocks);
}

Matrix make_H_full() {
    
    // Block diagonal: H_joint placed at rows [j*3, j*3+3), cols [j*12, j*12+12)
    Matrix H_full(N_MEAS, N_STATE);
    Matrix H_joint = make_H_joint();

    for (int j = 0; j < N_JOINTS; ++j)
        
        H_full.set_block(j * 3, j * 12, H_joint);

    return H_full;
}

Matrix make_R_full(double r_noise) {
    
    // r_noise = σ² (position variance in m²)
    Matrix R(N_MEAS, N_MEAS);
    for (int i = 0; i < N_MEAS; ++i)
        R(i, i) = r_noise;
    return R;
}

Matrix make_P0_full(double p_pos, double p_vel, double p_acc, double p_jerk) {
    
    // Initial covariance: diagonal per state component
    // Position: moderate uncertainty from first measurement
    // V/A/Jerk: large uncertainty (unknown)
    
    Matrix P0(N_STATE, N_STATE);

    for (int joint = 0; joint < N_JOINTS; ++joint) {
        
        int base = joint * N_JOINT;
        
        for (int axis = 0; axis < 3; ++axis) {
            
            int ab = base + axis * N_AXIS;
            P0(ab + 0, ab + 0) = p_pos;
            P0(ab + 1, ab + 1) = p_vel;
            P0(ab + 2, ab + 2) = p_acc;
            P0(ab + 3, ab + 3) = p_jerk;
        }
    }
    return P0;
}

Matrix make_x0(const std::vector<double>& z0) {
   
    // Initialise state from first measurement
    // Positions from z0, velocities/accelerations/jerks zero
    // z0 layout: [pelvis_x, pelvis_y, pelvis_z, L5_x, L5_y, L5_z, ...]
    
    if (static_cast<int>(z0.size()) != N_MEAS)
        
        throw std::invalid_argument("make_x0: z0 must have " + std::to_string(N_MEAS) + " elements");

    Matrix x0(N_STATE, 1);

    for (int j = 0; j < N_JOINTS; ++j) {
        
        int state_base = j * N_JOINT;
        int meas_base = j * 3;

        x0(state_base + 0, 0) = z0[meas_base + 0];  // px
        x0(state_base + 4, 0) = z0[meas_base + 1];  // py
        x0(state_base + 8, 0) = z0[meas_base + 2];  // pz
    }

    return x0;
}