// lkf.cpp
// Team Thunder — Kalman Filter Milestone 2
// LINEAR KALMAN FILTER — Full Body 3D Gait Estimation

#include "matrix_utils.h"
#include "csv_reader.h"
#include "kalman_matrices.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <stdexcept>

// LKFState carries filter state between timesteps
struct LKFState {
    
    Matrix x;   // state estimate (276×1)
    Matrix P;   // covariance (276×276)
};

// lkf_predicts time update (Eq. 58-59)
static LKFState lkf_predict(const LKFState& prev,
    
    const Matrix& F,
    const Matrix& Q) {
    LKFState pred;
    pred.x = F * prev.x;
    pred.P = F * prev.P * F.transpose() + Q;
    pred.P.symmetrize();
    return pred;
}

// lkf_update for measurement correction (Eq. 60-62, Joseph form)
static LKFState lkf_update(const LKFState& pred,
    
    const Matrix& H,
    const Matrix& R,
    const Matrix& z) {
    Matrix y = z - H * pred.x;
    Matrix Ht = H.transpose();
    Matrix PHt = pred.P * Ht;
    Matrix S = H * PHt + R;
    S.symmetrize();

    Matrix Kt = S.cholesky_solve(PHt.transpose());
    Matrix K = Kt.transpose();

    LKFState upd;
    upd.x = pred.x + K * y;

    Matrix I = Matrix::identity(N_STATE);
    Matrix IKH = I - K * H;
    upd.P = IKH * pred.P * IKH.transpose() + K * R * K.transpose();
    upd.P.symmetrize();

    return upd;
}

// CSV output helpers
static void write_csv_header(std::ofstream& out) {
    out << "frame";
    const char* axes[] = { "x", "y", "z" };
    const char* states[] = { "p", "v", "a", "j" };
    for (int j = 0; j < N_JOINTS; ++j)
        
        for (int ax = 0; ax < 3; ++ax)
            
            for (int s = 0; s < N_AXIS; ++s)
                out << "," << JOINT_NAMES[j] << "_" << states[s] << axes[ax];
    out << "\n";
}

static void write_state_row(std::ofstream& out, int frame, const Matrix& x) {
    
    out << frame;
    for (int i = 0; i < N_STATE; ++i)
        
        out << "," << std::setprecision(8) << std::fixed << x(i, 0);
   
    out << "\n";
}

// run_lkf 
static void run_lkf(const GaitDataset& noisy,
    
    const std::string& output_file,
    double dt,
    double phi,
    double r_noise) {

    std::cout << "\n[LKF] Building system matrices...\n";
    std::cout << "      dt      = " << dt << " s\n";
    std::cout << "      phi (Φ) = " << phi << "\n";
    std::cout << "      r_noise = " << r_noise << " m²\n";

    Matrix F = make_F_full(dt);
    Matrix Q = make_Q_full(dt, phi);
    Matrix H = make_H_full();
    Matrix R = make_R_full(r_noise);

    std::vector<double> z0_vec = noisy.get_measurement(0);
    Matrix z0(N_MEAS, 1);
    for (int i = 0; i < N_MEAS; ++i) z0(i, 0) = z0_vec[i];

    LKFState state;
    state.x = make_x0(z0_vec);
    state.P = make_P0_full();

    std::ofstream out(output_file);
    if (!out.is_open())
        
        throw std::runtime_error("[LKF] Cannot open output file: " + output_file);
    
    write_csv_header(out);
    write_state_row(out, 0, state.x);

    std::cout << "[LKF] Running filter over " << noisy.num_frames() << " frames...\n";

    auto t_start = std::chrono::high_resolution_clock::now();
    int N = noisy.num_frames();

    for (int k = 1; k < N; ++k) {
       
        LKFState predicted = lkf_predict(state, F, Q);

        std::vector<double> zk_vec = noisy.get_measurement(k);
        Matrix zk(N_MEAS, 1);
        for (int i = 0; i < N_MEAS; ++i) zk(i, 0) = zk_vec[i];

        state = lkf_update(predicted, H, R, zk);
        write_state_row(out, k, state.x);

        if (k % 500 == 0) {
            
            auto t_now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(t_now - t_start).count();
            std::cout << "  [LKF] frame " << k << " / " << N
                << "  (" << std::setprecision(1) << std::fixed
                << elapsed << " s elapsed)\n";
        }
    }

    out.close();

    auto t_end = std::chrono::high_resolution_clock::now();
    double total = std::chrono::duration<double>(t_end - t_start).count();

    std::cout << "\n[LKF] Done. " << N << " frames processed in "
        << std::setprecision(2) << total << " s\n";
    std::cout << "[LKF] Output written to: " << output_file << "\n";
}

// main
int main(int argc, char* argv[]) {
    
    std::cout << "=========================================\n";
    std::cout << "  Team Thunder — Linear Kalman Filter\n";
    std::cout << "  Milestone 2\n";
    std::cout << "=========================================\n";

    if (argc < 3) {
        
        std::cerr << "\nUsage: lkf <noisy_csv> <output_csv> [dt] [phi] [r_noise]\n";
        std::cerr << "  noisy_csv  : input CSV with 69 position columns\n";
        std::cerr << "  output_csv : output CSV for filtered state vectors\n";
        std::cerr << "  dt         : sampling interval in seconds (default 0.01)\n";
        std::cerr << "  phi        : process noise density Φ       (default 1.0)\n";
        std::cerr << "  r_noise    : measurement noise variance m² (default 0.001)\n";
        return 1;
    }

    std::string noisy_file = argv[1];
    std::string output_file = argv[2];
    
    double dt = (argc > 3) ? std::stod(argv[3]) : 0.01;
    double phi = (argc > 4) ? std::stod(argv[4]) : 1.0;
    double r_noise = (argc > 5) ? std::stod(argv[5]) : 0.001;

    try {
        
        GaitDataset noisy = load_gait_csv(noisy_file, dt);
        noisy.print_summary();
        run_lkf(noisy, output_file, dt, phi, r_noise);
    }
    catch (const std::exception& e) {
       
        std::cerr << "\n[LKF] ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}