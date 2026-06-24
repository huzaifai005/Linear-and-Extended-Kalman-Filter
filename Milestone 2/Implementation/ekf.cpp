#include "matrix_utils.h"
#include "csv_reader.h"
#include "kalman_matrices.h"
#include "arctan2_approx.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <stdexcept>

// EKFState
struct EKFState {
    Matrix x;   // 276×1
    Matrix P;   // 276×276
};

// a function that wraps angle to (-π, π]
static double wrap_angle(double a) {
    const double PI     = manual_pi();
    const double TWO_PI = 2.0 * PI;
    while (a >  PI) a -= TWO_PI;
    while (a < -PI) a += TWO_PI;
    return a;
}

// wrap angle components of the 69x1 measurement vector
static void wrap_innovation(Matrix& y) {
    for (int j = 0; j < N_JOINTS; ++j) {
        y(j*3 + 1, 0) = wrap_angle(y(j*3 + 1, 0));   // azimuth
        y(j*3 + 2, 0) = wrap_angle(y(j*3 + 2, 0));   // elevation
    }
}

// Builds the 69×69 measurement noise covariance matrix per joint, per timestep.
static Matrix make_R_adaptive(const Matrix& x_pred, double sigma_cart) {
    const double sigma2   = sigma_cart * sigma_cart;
    const double min_rho2 = 0.01;   
    const double min_r2   = 0.01;

    Matrix R(N_MEAS, N_MEAS);  

    for (int j = 0; j < N_JOINTS; ++j) {
        double px = x_pred(j*N_JOINT + 0, 0);
        double py = x_pred(j*N_JOINT + 4, 0);
        double pz = x_pred(j*N_JOINT + 8, 0);

        double rho2 = px*px + py*py;
        double r2   = rho2 + pz*pz;

        if (rho2 < min_rho2) rho2 = min_rho2;
        if (r2   < min_r2  ) r2   = min_r2;

        int base = j * 3;
        R(base + 0, base + 0) = sigma2;           // range variance
        R(base + 1, base + 1) = sigma2 / rho2;    // azimuth variance
        R(base + 2, base + 2) = sigma2 / r2;      // elevation variance
    }
    return R;
}

// conversion of cartesian positions to spherical ones
static Matrix cartesian_to_spherical(const std::vector<double>& z_cart) {
    Matrix z_sph(N_MEAS, 1);
    for (int j = 0; j < N_JOINTS; ++j) {
        int base = j * 3;
        double px = z_cart[base + 0];
        double py = z_cart[base + 1];
        double pz = z_cart[base + 2];
        double rho = manual_sqrt(px*px + py*py);
        double r   = manual_sqrt(px*px + py*py + pz*pz);
        static const double EPS = 1e-9;
        z_sph(base + 0, 0) = r;
        z_sph(base + 1, 0) = (r > EPS) ? manual_atan2(py, px) : 0.0;
        z_sph(base + 2, 0) = (r > EPS) ? manual_atan2(pz, rho) : 0.0;
    }
    return z_sph;
}

// time update
static EKFState ekf_predict(const EKFState& prev,
                              const Matrix& F,
                              const Matrix& Q) {
    EKFState pred;
    pred.x = F * prev.x;
    pred.P = F * prev.P * F.transpose() + Q;
    pred.P.symmetrize();
    return pred;
}

// update step
static EKFState ekf_update(const EKFState& pred,
                             double sigma_cart,
                             const Matrix& z_sph) {

    Matrix R = make_R_adaptive(pred.x, sigma_cart);

    // Jacobian at predicted state
    Matrix Hk  = H_jacobian_full(pred.x); 
    Matrix Hkt = Hk.transpose();             

    // Predicted measurement
    Matrix z_hat = h_full(pred.x);          

    // Innovation with angle wrapping
    Matrix y = z_sph - z_hat;
    wrap_innovation(y);

    
    Matrix PHkt = pred.P * Hkt;
    Matrix S    = Hk * PHkt + R;
    S.symmetrize();

    // Kalman gain via Cholesky solve
    Matrix Kt = S.cholesky_solve(PHkt.transpose());
    Matrix K  = Kt.transpose();     

    // State update
    EKFState upd;
    upd.x = pred.x + K * y;

    // Covariance - Joseph form
    Matrix I   = Matrix::identity(N_STATE);
    Matrix IKH = I - K * Hk;
    upd.P = IKH * pred.P * IKH.transpose()
          + K   * R      * K.transpose();
    upd.P.symmetrize();

    return upd;
}

// writes output file with column names per joint and axis
// and state vector rows frame by frame.
static void write_csv_header(std::ofstream& out) {
    out << "frame";
    const char* axes[]   = {"x","y","z"};
    const char* states[] = {"p","v","a","j"};
    for (int j = 0; j < N_JOINTS; ++j)
        for (int ax = 0; ax < 3; ++ax)
            for (int s = 0; s < N_AXIS; ++s)
                out << "," << JOINT_NAMES[j]
                    << "_" << states[s] << axes[ax];
    out << "\n";
}

static void write_state_row(std::ofstream& out, int frame, const Matrix& x) {
    out << frame;
    for (int i = 0; i < N_STATE; ++i)
        out << "," << std::setprecision(8) << std::fixed << x(i,0);
    out << "\n";
}

// run_ekf
static void run_ekf(const GaitDataset& noisy,
                    const std::string& output_file,
                    double dt,
                    double phi,
                    double sigma_cart) {

    std::cout << "\n[EKF] Building system matrices...\n";
    std::cout << "      dt         = " << dt         << " s\n";
    std::cout << "      phi (Φ)    = " << phi        << "\n";
    std::cout << "      sigma_cart = " << sigma_cart << " m\n";
    std::cout << "      R          = adaptive (computed per step per joint)\n";

    Matrix F = make_F_full(dt);
    Matrix Q = make_Q_full(dt, phi);

    std::cout << "[EKF] F: 276×276   Q: 276×276   Hk: 69×276 (per step)\n";

    // Initialise
    std::vector<double> z0_vec = noisy.get_measurement(0);
    EKFState state;
    state.x = make_x0(z0_vec);
    state.P = make_P0_full();

    // Output file
    std::ofstream out(output_file);
    if (!out.is_open())
        throw std::runtime_error("[EKF] Cannot open output: " + output_file);
    write_csv_header(out);
    write_state_row(out, 0, state.x);

    std::cout << "[EKF] Running filter over "
              << noisy.num_frames() << " frames...\n";

    auto t_start = std::chrono::high_resolution_clock::now();
    int N = noisy.num_frames();

    for (int k = 1; k < N; ++k) {

        // Predict
        EKFState predicted = ekf_predict(state, F, Q);

        // Convert Cartesian measurement to spherical
        std::vector<double> zk_cart = noisy.get_measurement(k);
        Matrix zk_sph = cartesian_to_spherical(zk_cart);

        // Update
        state = ekf_update(predicted, sigma_cart, zk_sph);

        write_state_row(out, k, state.x);

        if (k % 500 == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            double el = std::chrono::duration<double>(now - t_start).count();
            std::cout << "  [EKF] frame " << k << "/" << N
                      << "  (" << std::fixed << std::setprecision(1)
                      << el << "s)\n";
        }
    }

    out.close();

    double total = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    std::cout << "\n[EKF] Done. " << N << " frames in "
              << std::setprecision(2) << total << "s\n";
    std::cout << "[EKF] Output: " << output_file << "\n";
}

// main function to test
int main(int argc, char* argv[]) {
    std::cout << "=========================================\n";
    std::cout << "  Team Thunder — Extended Kalman Filter\n";
    std::cout << "  Milestone 2\n";
    std::cout << "=========================================\n";

    if (argc < 3) {
        std::cerr << "\nUsage: ekf <noisy_csv> <output_csv> [dt] [phi] [sigma_cart]\n";
        std::cerr << "  dt         : sampling interval (default 0.01 s)\n";
        std::cerr << "  phi        : process noise Φ   (default 0.01)\n";
        std::cerr << "  sigma_cart : Cartesian noise std dev in m (default 0.05)\n";
        std::cerr << "\nExample:\n";
        std::cerr << "  ./ekf noisy_data.csv ekf_output.csv 0.01 0.01 0.05\n\n";
        return 1;
    }

    std::string noisy_file  = argv[1];
    std::string output_file = argv[2];
    double dt         = (argc > 3) ? std::stod(argv[3]) : 0.01;
    double phi        = (argc > 4) ? std::stod(argv[4]) : 0.01;
    double sigma_cart = (argc > 5) ? std::stod(argv[5]) : 0.05;

    try {
        GaitDataset noisy = load_gait_csv(noisy_file, dt);
        noisy.print_summary();
        run_ekf(noisy, output_file, dt, phi, sigma_cart);
    } catch (const std::exception& e) {
        std::cerr << "\n[EKF] ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
