#include "matrix_utils.h"
#include "csv_reader.h"
#include "kalman_matrices.h"
#include "arctan2_approx.h"
#include <iostream>
#include <cmath> 
#include <cassert>
#include <vector>

// it is used to verify the numeric correctness
static bool approx_equal(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) < tol;
}

// test functions to validate different operations
void test_basic_arithmetic() {
    std::cout << "---- Test 1: Basic arithmetic ----\n";
    Matrix A(2,2); A(0,0)=1;A(0,1)=2;A(1,0)=3;A(1,1)=4;
    Matrix B(2,2); B(0,0)=5;B(0,1)=6;B(1,0)=7;B(1,1)=8;
    Matrix D = A * B;
    assert(approx_equal(D(0,0),19)); assert(approx_equal(D(1,1),50));
    std::cout << "  PASS\n";
}
void test_transpose() {
    std::cout << "---- Test 2: Transpose ----\n";
    Matrix A(2,3); A(0,0)=1;A(0,1)=2;A(0,2)=3;A(1,0)=4;A(1,1)=5;A(1,2)=6;
    Matrix At = A.transpose();
    assert(At.rows==3 && At.cols==2);
    assert(approx_equal(At(2,1),6));
    std::cout << "  PASS\n";
}
void test_cholesky_solve() {
    std::cout << "---- Test 3: Cholesky solve ----\n";
    Matrix A(3,3);
    A(0,0)=4;A(0,1)=2;A(0,2)=2;
    A(1,0)=2;A(1,1)=3;A(1,2)=1;
    A(2,0)=2;A(2,1)=1;A(2,2)=3;
    Matrix B(3,2); B(0,0)=1;B(1,1)=1;
    Matrix X = A.cholesky_solve(B);
    Matrix AX = A * X;
    for(int i=0;i<3;++i) for(int j=0;j<2;++j)
        assert(approx_equal(AX(i,j),B(i,j),1e-9));
    std::cout << "  PASS\n";
}
void test_F1_kinematics() {
    std::cout << "---- Test 4: F1 kinematics ----\n";
    double dt=0.1; Matrix F1=make_F1(dt);
    Matrix x(4,1); x(0,0)=0;x(1,0)=1;x(2,0)=0;x(3,0)=0;
    Matrix xn=F1*x;
    assert(approx_equal(xn(0,0),0.1,1e-10));
    assert(approx_equal(xn(1,0),1.0,1e-10));
    std::cout << "  PASS\n";
}
void test_Q1_properties() {
    std::cout << "---- Test 5: Q1 SPD ----\n";
    Matrix Q1=make_Q1(0.01,1.0);
    assert(Q1.is_symmetric(1e-12));
    Matrix L=Q1.cholesky(); (void)L;
    std::cout << "  PASS\n";
}
void test_fullbody_dimensions() {
    std::cout << "---- Test 6: Full-body matrix dimensions ----\n";
    Matrix F=make_F_full(0.01);
    Matrix Q=make_Q_full(0.01,1.0);
    Matrix H=make_H_full();
    assert(F.rows==276&&F.cols==276);
    assert(Q.rows==276&&Q.cols==276);
    assert(H.rows== 69&&H.cols==276);
    std::cout << "  F:276×276  Q:276×276  H:69×276  PASS\n";
}
void test_H_extracts_positions() {
    std::cout << "---- Test 7: H extracts positions ----\n";
    Matrix H=make_H_full();
    Matrix x(276,1);
    x(0,0)=1.5;x(4,0)=0.8;x(8,0)=0.9;
    Matrix z=H*x;
    assert(approx_equal(z(0,0),1.5));
    assert(approx_equal(z(1,0),0.8));
    assert(approx_equal(z(2,0),0.9));
    std::cout << "  PASS\n";
}
void test_single_lkf_step() {
    std::cout << "---- Test 8: Single LKF step ----\n";
    Matrix F=make_F_full(0.01),Q=make_Q_full(0.01,1.0);
    Matrix H=make_H_full(),R=make_R_full(0.001);
    Matrix I=Matrix::identity(N_STATE);
    std::vector<double> z0(N_MEAS,0.0); z0[0]=1.0;z0[1]=0.5;z0[2]=0.8;
    Matrix x=make_x0(z0); Matrix P=make_P0_full();
    Matrix x_pred=F*x; Matrix P_pred=F*P*F.transpose()+Q; P_pred.symmetrize();
    std::vector<double> z1v(N_MEAS,0.0); z1v[0]=1.01;z1v[1]=0.5;z1v[2]=0.8;
    Matrix z1(N_MEAS,1); for(int i=0;i<N_MEAS;++i) z1(i,0)=z1v[i];
    Matrix Ht=H.transpose(),PHt=P_pred*Ht,S=H*PHt+R; S.symmetrize();
    Matrix K=S.cholesky_solve(PHt.transpose()).transpose();
    Matrix IKH=I-K*H;
    Matrix P_upd=IKH*P_pred*IKH.transpose()+K*R*K.transpose(); P_upd.symmetrize();
    assert(P_upd(0,0) < P_pred(0,0));
    assert(P_upd.is_symmetric(1e-8));
    std::cout << "  PASS\n";
}

// this checks if our manual implementation is actually correct or not
void test_manual_sqrt() {
    std::cout << "---- Test 9: manual_sqrt accuracy ----\n";
    double cases[] = {0.0, 0.001, 0.25, 1.0, 2.0, 4.0, 100.0, 12345.6789};
    for (double s : cases) {
        double ref = std::sqrt(s);
        double our = manual_sqrt(s);
        double err = std::abs(our - ref);
        double tol = (ref > 1e-6) ? 1e-12 * ref : 1e-12;
        if (err > tol) {
            std::cerr << "  FAIL: sqrt(" << s << ") ref=" << ref
                      << " our=" << our << " err=" << err << "\n";
            assert(false);
        }
    }
    std::cout << "  PASS: manual_sqrt matches std::sqrt to 1e-12 relative error\n";
}

void test_manual_atan2() {
    std::cout << "---- Test 10: manual_atan2 accuracy (CORDIC) ----\n";

    // {y, x, description}
    struct Case { double y, x; const char* desc; };
    Case cases[] = {
        { 0.0,   1.0,  "positive x-axis (0)" },
        { 1.0,   1.0,  "Q1 45°" },
        { 1.0,   0.0,  "positive y-axis (π/2)" },
        { 1.0,  -1.0,  "Q2 135°" },
        { 0.0,  -1.0,  "negative x-axis (π)" },
        {-1.0,  -1.0,  "Q3 -135°" },
        {-1.0,   0.0,  "negative y-axis (-π/2)" },
        {-1.0,   1.0,  "Q4 -45°" },
        { 3.0,   4.0,  "Q1 non-trivial" },
        {-3.0,   4.0,  "Q4 non-trivial" },
        { 3.0,  -4.0,  "Q2 non-trivial" },
        {-3.0,  -4.0,  "Q3 non-trivial" },
        { 0.001, 1.0,  "near x-axis" },
        { 1.0,  0.001, "near y-axis" },
    };

    double max_err = 0.0;
    const char* worst = "";

    for (auto& c : cases) {
        double ref = std::atan2(c.y, c.x);
        double our = manual_atan2(c.y, c.x);
        double err = std::abs(our - ref);
        if (err > max_err) { max_err = err; worst = c.desc; }
        if (err > 1e-6) {
            std::cerr << "  FAIL at (" << c.desc << "): "
                      << "ref=" << ref << " our=" << our
                      << " err=" << err << "\n";
            assert(false);
        }
    }

    std::cout << "  PASS: max error = " << max_err
              << " (at '" << worst << "')\n";
    std::cout << "  CORDIC (" << CORDIC_ITERATIONS
              << " iterations) matches std::atan2 to ~1e-" 
              << (int)(-std::log10(max_err + 1e-20))
              << " precision\n";
}

// check known values
void test_h_joint() {
    std::cout << "---- Test 11: h_joint nonlinear measurement ----\n";

    Matrix x_j(12, 1);
    x_j(0,0) = 3.0;   // px
    x_j(4,0) = 4.0;   // py
    x_j(8,0) = 0.0;   // pz

    Matrix z = h_joint(x_j);

    double r_ref     = 5.0;
    double theta_ref = std::atan2(4.0, 3.0);
    double phi_ref   = 0.0;

    assert(approx_equal(z(0,0), r_ref,     1e-9));
    assert(approx_equal(z(1,0), theta_ref, 1e-6));
    assert(approx_equal(z(2,0), phi_ref,   1e-6));

    std::cout << "  r=" << z(0,0) << " θ=" << z(1,0)
              << " φ=" << z(2,0) << "\n";

    Matrix x_j2(12,1);
    x_j2(8,0) = 5.0;   // pz only
    Matrix z2 = h_joint(x_j2);
    assert(approx_equal(z2(0,0), 5.0,   1e-9));

    std::cout << "  PASS\n";
}

// very important to verify otherwise the filter would diverge
void test_jacobian_numerical() {
    std::cout << "---- Test 12: Jacobian numerical verification ----\n";

    // Generic non-degenerate position
    Matrix x_j(12,1);
    x_j(0,0) = 1.2;   // px
    x_j(4,0) = 0.5;   // py
    x_j(8,0) = 0.8;   // pz

    Matrix H_analytic = H_jacobian_joint(x_j);

    auto h_ref = [](const Matrix& xj) -> Matrix {
        double px = xj(0,0), py = xj(4,0), pz = xj(8,0);
        double rho = std::sqrt(px*px + py*py);
        double r   = std::sqrt(px*px + py*py + pz*pz);
        Matrix z(3,1);
        z(0,0) = r;
        z(1,0) = std::atan2(py, px);
        z(2,0) = std::atan2(pz, rho);
        return z;
    };

    // Central-difference numerical Jacobian
    double eps = 1e-5;
    Matrix H_numerical(3, 12);
    for (int col = 0; col < 12; ++col) {
        Matrix xp = x_j, xm = x_j;
        xp(col,0) += eps;
        xm(col,0) -= eps;
        Matrix zp = h_ref(xp), zm = h_ref(xm);
        for (int row = 0; row < 3; ++row)
            H_numerical(row, col) = (zp(row,0) - zm(row,0)) / (2.0*eps);
    }

    // Comparison
    double max_err = 0.0;
    int worst_r = 0, worst_c = 0;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 12; ++c) {
            double err = std::abs(H_analytic(r,c) - H_numerical(r,c));
            if (err > max_err) { max_err = err; worst_r=r; worst_c=c; }
        }

    std::cout << "  Max |analytic - numerical| = " << max_err
              << "  (at row=" << worst_r << " col=" << worst_c << ")\n";

    if (max_err > 1e-6) {
        std::cerr << "  FAIL: Jacobian mismatch exceeds tolerance\n";
        H_analytic.print("H_analytic", 3);
        H_numerical.print("H_numerical", 3);
        assert(false);
    }
    std::cout << "  PASS: analytic Jacobian matches finite-difference to 1e-6\n";
}

// here we test the full body measurement model
void test_full_ekf_dimensions() {
    std::cout << "---- Test 13: h_full and H_jacobian_full dimensions ----\n";

    Matrix x(276,1);
    // Set non-zero positions to avoid degenerate singularities
    for (int j = 0; j < N_JOINTS; ++j) {
        x(j*12+0, 0) = 0.5 + j*0.1;   // px
        x(j*12+4, 0) = 0.3 + j*0.05;  // py
        x(j*12+8, 0) = 0.8 + j*0.02;  // pz
    }

    Matrix z  = h_full(x);
    Matrix Hk = H_jacobian_full(x);

    assert(z.rows  == 69  && z.cols  == 1);
    assert(Hk.rows == 69  && Hk.cols == 276);

    for (int r = 3; r < 6; ++r)
        for (int c = 0; c < 12; ++c)
            assert(approx_equal(Hk(r,c), 0.0, 1e-15));

    std::cout << "  z:  69x1   Hk: 69x276   block-diagonal structure ✓\n";
    std::cout << "  PASS\n";
}

// here we simulate one complete ekf cycle
void test_single_ekf_step() {
    std::cout << "---- Test 14: Single EKF predict+update step ----\n";

    double dt=0.01, phi=1.0, r_range=0.001, r_angle=0.01;

    Matrix F = make_F_full(dt);
    Matrix Q = make_Q_full(dt, phi);
    Matrix I = Matrix::identity(N_STATE);

    // R for EKF that block diagonal diag(r_range, r_angle, r_angle) × 23
    Matrix R(N_MEAS, N_MEAS);
    for (int j = 0; j < N_JOINTS; ++j) {
        R(j*3+0,j*3+0) = r_range;
        R(j*3+1,j*3+1) = r_angle;
        R(j*3+2,j*3+2) = r_angle;
    }

    // Initialise that pelvis at (1.0, 0.5, 0.8), all others at small offsets
    std::vector<double> z0(N_MEAS, 0.0);
    for (int j = 0; j < N_JOINTS; ++j) {
        z0[j*3+0] = 0.5 + j*0.1;
        z0[j*3+1] = 0.3;
        z0[j*3+2] = 0.8;
    }
    Matrix x = make_x0(z0);
    Matrix P = make_P0_full();

    // Predict
    Matrix x_pred = F * x;
    Matrix P_pred = F * P * F.transpose() + Q;
    P_pred.symmetrize();

    // Linearise
    Matrix Hk  = H_jacobian_full(x_pred);  
    Matrix Hkt = Hk.transpose();
    Matrix z_hat = h_full(x_pred);          

    // Measurement: same as z0 converted to spherical
    Matrix z_meas(N_MEAS, 1);
    for (int j = 0; j < N_JOINTS; ++j) {
        double px = z0[j*3+0], py = z0[j*3+1], pz = z0[j*3+2];
        double rho = manual_sqrt(px*px + py*py);
        double r   = manual_sqrt(px*px + py*py + pz*pz);
        z_meas(j*3+0, 0) = r;
        z_meas(j*3+1, 0) = manual_atan2(py, px);
        z_meas(j*3+2, 0) = manual_atan2(pz, rho);
    }

    // Innovation
    Matrix y = z_meas - z_hat;

    // S, K
    Matrix PHkt = P_pred * Hkt;
    Matrix S    = Hk * PHkt + R;
    S.symmetrize();
    Matrix K = S.cholesky_solve(PHkt.transpose()).transpose();

    // Update
    Matrix x_upd = x_pred + K * y;
    Matrix IKH   = I - K * Hk;
    Matrix P_upd = IKH * P_pred * IKH.transpose()
                 + K   * R      * K.transpose();
    P_upd.symmetrize();

    // Checks
    double p_before = P_pred(0,0);
    double p_after  = P_upd(0,0);
    std::cout << "  P[0,0]: before=" << p_before << "  after=" << p_after << "\n";
    assert(p_after < p_before);               // covariance must decrease
    assert(P_upd.is_symmetric(1e-8));          // P must stay symmetric
    assert(x_upd.rows == 276);                 // correct dimensions

    std::cout << "  PASS: EKF step — Jacobian applied, covariance reduced, P symmetric\n";
}

void test_csv(const std::string& tf, const std::string& nf) {
    std::cout << "---- Test 15: CSV loading ----\n";
    auto [td, nd] = load_both_datasets(tf, nf);
    td.print_summary(); nd.print_summary();
    assert(static_cast<int>(td.get_measurement(0).size()) == N_MEAS);
    std::cout << "  PASS\n";
}

// main
int main(int argc, char* argv[]) {
    std::cout << "=========================================\n";
    std::cout << "  Team Thunder — Milestone 2 Unit Tests\n";
    std::cout << "  Phase 3+4+5: Matrices + LKF + EKF\n";
    std::cout << "=========================================\n\n";

    try {
        // Matrix utils
        test_basic_arithmetic();
        test_transpose();
        test_cholesky_solve();

        // Kalman matrices
        test_F1_kinematics();
        test_Q1_properties();
        test_fullbody_dimensions();
        test_H_extracts_positions();
        test_single_lkf_step();

        // arctan2 / EKF
        test_manual_sqrt();
        test_manual_atan2();
        test_h_joint();
        test_jacobian_numerical();
        test_full_ekf_dimensions();
        test_single_ekf_step();

        std::cout << "\n✓  All 14 tests PASSED.\n\n";

        if (argc == 3) {
            test_csv(argv[1], argv[2]);
            std::cout << "\n✓  CSV test PASSED.\n";
        } else {
            std::cout << "  (No CSV files provided — run with:\n";
            std::cout << "   ./test_main true_data.csv noisy_data.csv)\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "\n✗  TEST FAILED: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n  Both LKF and EKF are ready to run.\n";
    return 0;
}
