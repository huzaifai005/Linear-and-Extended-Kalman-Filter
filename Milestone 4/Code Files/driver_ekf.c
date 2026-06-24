#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N_JOINTS 23
#define N_AXIS 4
#define N_JOINT 12
#define N_STATE 276
#define N_MEAS 69
#define MAX_FRAMES 4000

#define DT 0.01
#define PHI 0.01
#define SIGMA_CART 0.707

#define M(A, i, j, cols) ((A)[(long)(i) * (cols) + (j)])

static const char *JOINT_NAMES[N_JOINTS] = {
    "pelvis", "L5", "L3", "T12", "T8", "neck", "head",
    "shoulderRight", "upperArmRight", "forearmRight", "handRight",
    "shoulderLeft", "upperArmLeft", "forearmLeft", "handLeft",
    "upperLegRight", "lowerLegRight", "footRight", "toeRight",
    "upperLegLeft", "lowerLegLeft", "footLeft", "toeLeft"};

// Assembly declarations
extern void asm_mat_mul(double *, double *, double *, long, long, long);
extern void asm_mat_add(double *, double *, double *, long, long);
extern void asm_mat_sub(double *, double *, double *, long, long);
extern void asm_mat_transpose(double *, double *, long, long);
extern void asm_symmetrize(double *, long);
extern void asm_cholesky_solve(double *, double *, double *, double *, long, long);
extern double asm_newton_sqrt(double);
extern void asm_h_full(double *, double *);
extern void asm_jacobian_full(double *, double *, double *);
extern void asm_vec_zero_matrix(double *, long);
extern void asm_vec_mat_copy(double *, double *, long);
extern double asm_cordic_atan2(double, double);
extern double asm_wrap_angle(double);

// Helpers
static double *alloc_mat(long r, long c)
{
    double *p = (double *)calloc(r * c, sizeof(double));
    if (!p)
    {
        fprintf(stderr, "OOM\n");
        exit(1);
    }
    return p;
}
static void mat_copy(double *d, double *s, long n) { memcpy(d, s, n * sizeof(double)); }
static void mat_negate(double *A, long n)
{
    for (long i = 0; i < n; i++)
        A[i] = -A[i];
}
static void add_identity(double *A, long n)
{
    for (long i = 0; i < n; i++)
        M(A, i, i, n) += 1.0;
}

// Build 4×4 F1 block
static void build_F1(double *F1)
{
    double dt = DT, dt2 = dt * dt, dt3 = dt2 * dt;
    memset(F1, 0, 16 * sizeof(double));
    M(F1, 0, 0, 4) = 1;
    M(F1, 0, 1, 4) = dt;
    M(F1, 0, 2, 4) = dt2 / 2;
    M(F1, 0, 3, 4) = dt3 / 6;
    M(F1, 1, 1, 4) = 1;
    M(F1, 1, 2, 4) = dt;
    M(F1, 1, 3, 4) = dt2 / 2;
    M(F1, 2, 2, 4) = 1;
    M(F1, 2, 3, 4) = dt;
    M(F1, 3, 3, 4) = 1;
}
static void build_Q1(double *Q1, double phi)
{
    double dt = DT, dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt3 * dt,
           dt5 = dt4 * dt, dt6 = dt5 * dt, dt7 = dt6 * dt;
    memset(Q1, 0, 16 * sizeof(double));
    M(Q1, 0, 0, 4) = phi * dt7 / 252;
    M(Q1, 0, 1, 4) = phi * dt6 / 72;
    M(Q1, 1, 0, 4) = M(Q1, 0, 1, 4);
    M(Q1, 0, 2, 4) = phi * dt5 / 30;
    M(Q1, 2, 0, 4) = M(Q1, 0, 2, 4);
    M(Q1, 0, 3, 4) = phi * dt4 / 24;
    M(Q1, 3, 0, 4) = M(Q1, 0, 3, 4);
    M(Q1, 1, 1, 4) = phi * dt5 / 20;
    M(Q1, 1, 2, 4) = phi * dt4 / 8;
    M(Q1, 2, 1, 4) = M(Q1, 1, 2, 4);
    M(Q1, 1, 3, 4) = phi * dt3 / 6;
    M(Q1, 3, 1, 4) = M(Q1, 1, 3, 4);
    M(Q1, 2, 2, 4) = phi * dt3 / 3;
    M(Q1, 2, 3, 4) = phi * dt2 / 2;
    M(Q1, 3, 2, 4) = M(Q1, 2, 3, 4);
    M(Q1, 3, 3, 4) = phi * dt;
}
static void build_P0(double *P)
{
    memset(P, 0, (long)N_STATE * N_STATE * sizeof(double));
    for (int j = 0; j < N_JOINTS; j++)
    {
        int b = j * N_JOINT;
        for (int a = 0; a < 3; a++)
        {
            int ab = b + a * N_AXIS;
            M(P, ab + 0, ab + 0, N_STATE) = 1.0;
            M(P, ab + 1, ab + 1, N_STATE) = 100.0;
            M(P, ab + 2, ab + 2, N_STATE) = 100.0;
            M(P, ab + 3, ab + 3, N_STATE) = 100.0;
        }
    }
}
static void build_x0(double *x, double *z0)
{
    memset(x, 0, N_STATE * sizeof(double));
    for (int j = 0; j < N_JOINTS; j++)
    {
        x[j * N_JOINT + 0] = z0[j * 3 + 0];
        x[j * N_JOINT + 4] = z0[j * 3 + 1];
        x[j * N_JOINT + 8] = z0[j * 3 + 2];
    }
}

// Block-diagonal predict
static void ekf_predict_fast(double *x, double *P,
                             double *F1, double *Q1, double *F1T,
                             double *x_pred, double *P_pred)
{
    for (int b = 0; b < 69; b++)
    {
        int row = b * 4;
        double xb[4] = {x[row], x[row + 1], x[row + 2], x[row + 3]};
        double yb[4];
        asm_mat_mul(F1, xb, yb, 4, 4, 1);
        x_pred[row] = yb[0];
        x_pred[row + 1] = yb[1];
        x_pred[row + 2] = yb[2];
        x_pred[row + 3] = yb[3];
    }

    memset(P_pred, 0, (long)N_STATE * N_STATE * sizeof(double));
    double tmp4[16], fp[16], fpft[16];
    for (int bi = 0; bi < 69; bi++)
    {
        for (int bj = 0; bj < 69; bj++)
        {
            int ri = bi * 4, rj = bj * 4;
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < 4; c++)
                    tmp4[r * 4 + c] = M(P, ri + r, rj + c, N_STATE);
            asm_mat_mul(F1, tmp4, fp, 4, 4, 4);
            asm_mat_mul(fp, F1T, fpft, 4, 4, 4);
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < 4; c++)
                    M(P_pred, ri + r, rj + c, N_STATE) += fpft[r * 4 + c];
        }
    }
    for (int b = 0; b < 69; b++)
    {
        int base = b * 4;
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                M(P_pred, base + r, base + c, N_STATE) += Q1[r * 4 + c];
    }
    asm_symmetrize(P_pred, N_STATE);
}

// Cartesian -> spherical
static void cartesian_to_spherical(double *z_cart, double *z_sph)
{
    for (int j = 0; j < N_JOINTS; j++)
    {
        double px = z_cart[j * 3 + 0], py = z_cart[j * 3 + 1], pz = z_cart[j * 3 + 2];
        double rho = sqrt(px * px + py * py);
        double r   = sqrt(px * px + py * py + pz * pz);
        double eps = 1e-9;
        z_sph[j * 3 + 0] = r;
        z_sph[j * 3 + 1] = (r > eps) ? atan2(py, px) : 0.0;
        z_sph[j * 3 + 2] = (r > eps) ? atan2(pz, rho) : 0.0;
    }
}

// Adaptive R (diagonal, recomputed each frame from x_pred)
static void build_R_adaptive(double *R, double *x_pred, double sigma_cart)
{
    double sigma2 = sigma_cart * sigma_cart, min_rho2 = 0.01;
    memset(R, 0, (long)N_MEAS * N_MEAS * sizeof(double));
    for (int j = 0; j < N_JOINTS; j++)
    {
        double px = x_pred[j * N_JOINT + 0], py = x_pred[j * N_JOINT + 4], pz = x_pred[j * N_JOINT + 8];
        double rho2 = px * px + py * py;
        if (rho2 < min_rho2)
            rho2 = min_rho2;
        double r2 = rho2 + pz * pz;
        if (r2 < min_rho2)
            r2 = min_rho2;
        int b = j * 3;
        M(R, b + 0, b + 0, N_MEAS) = sigma2;
        M(R, b + 1, b + 1, N_MEAS) = sigma2 / rho2;
        M(R, b + 2, b + 2, N_MEAS) = sigma2 / r2;
    }
}

// EKF update
static void ekf_update(double *x_pred, double *P_pred,
                       double *Hk, double *R, double *z_sph,
                       double *x_upd, double *P_upd,
                       double *S, double *PHt, double *Kt,
                       double *z_hat, double *y,
                       double *IKH, double *tmp1, double *tmp2,
                       double *Ychol, double *K)
{
    // PHt = P_pred * Hk^T  (276x276 x 276x69 = 276x69)
    asm_mat_transpose(Hk, IKH, N_MEAS, N_STATE);
    asm_mat_mul(P_pred, IKH, PHt, N_STATE, N_STATE, N_MEAS);

    // S = Hk * PHt + R  (69x276 x 276x69 = 69x69)
    asm_mat_mul(Hk, PHt, S, N_MEAS, N_STATE, N_MEAS);
    asm_mat_add(S, R, S, N_MEAS, N_MEAS);
    asm_symmetrize(S, N_MEAS);

    // Cholesky solve: S * K^T = PHt^T
    asm_mat_transpose(PHt, Kt, N_STATE, N_MEAS);
    asm_vec_mat_copy(Ychol, S, (long)N_MEAS * N_MEAS);
    asm_cholesky_solve(Ychol, Kt, Kt, tmp1, N_MEAS, N_STATE);
    asm_mat_transpose(Kt, K, N_MEAS, N_STATE);

    // y = z_sph - h_full(x_pred), angle-wrapped
    asm_h_full(x_pred, z_hat);
    asm_mat_sub(z_sph, z_hat, y, N_MEAS, 1);
    for (int j = 0; j < N_JOINTS; j++)
    {
        y[j * 3 + 1] = asm_wrap_angle(y[j * 3 + 1]);
        y[j * 3 + 2] = asm_wrap_angle(y[j * 3 + 2]);
    }

    // x_upd = x_pred + K*y
    asm_mat_mul(K, y, x_upd, N_STATE, N_MEAS, 1);
    asm_mat_add(x_pred, x_upd, x_upd, N_STATE, 1);

    // P_upd = P_pred - K*(Hk*P_pred)

    asm_mat_mul(Hk,  P_pred, tmp1, N_MEAS,  N_STATE, N_STATE);
    asm_mat_mul(K,   tmp1,   tmp2, N_STATE, N_MEAS,  N_STATE);
    asm_mat_sub(P_pred, tmp2, P_upd, N_STATE, N_STATE);
    asm_symmetrize(P_upd, N_STATE);
}

// CSV Load
static int load_csv(const char *fn, double data[][N_MEAS], int *nf)
{
    FILE *f = fopen(fn, "r");
    if (!f)
    {
        perror(fn);
        return -1;
    }
    char line[8192];
    int frame = 0, hdr = 0;
    while (fgets(line, sizeof(line), f))
    {
        if (!hdr)
        {
            char *p = line;
            while (*p == ' ' || *p == '\t')
                p++;
            if (*p < '0' || *p > '9')
            {
                hdr = 1;
                continue;
            }
            hdr = 1;
        }
        if (frame >= MAX_FRAMES)
            break;
        char *tok = strtok(line, ",\n\r");
        for (int i = 0; i < N_MEAS && tok; i++)
        {
            data[frame][i] = atof(tok);
            tok = strtok(NULL, ",\n\r");
        }
        frame++;
    }
    fclose(f);
    *nf = frame;
    printf("[driver_ekf] Loaded %d frames from '%s'\n", frame, fn);
    return 0;
}

// CSV write 
static void write_output(const char *fn, double **states, int nf)
{
    FILE *f = fopen(fn, "w");
    if (!f)
    {
        perror(fn);
        return;
    }
    fprintf(f, "frame");
    const char *ax[] = {"x", "y", "z"}, *st[] = {"p", "v", "a", "j"};
    for (int j = 0; j < N_JOINTS; j++)
        for (int a = 0; a < 3; a++)
            for (int s = 0; s < N_AXIS; s++)
                fprintf(f, ",%s_%s%s", JOINT_NAMES[j], st[s], ax[a]);
    fprintf(f, "\n");
    for (int k = 0; k < nf; k++)
    {
        fprintf(f, "%d", k);
        for (int i = 0; i < N_STATE; i++)
            fprintf(f, ",%.8f", states[k][i]);
        fprintf(f, "\n");
    }
    fclose(f);
    printf("[driver_ekf] Output written to '%s'\n", fn);
}

// main
int main(int argc, char *argv[])
{
    printf("=========================================\n");
    printf("  Team Thunder — EKF Assembly Driver\n");
    printf("  Milestone 3 — Block-Diagonal Optimised\n");
    printf("=========================================\n\n");
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s <noisy.csv> <out.csv>\n", argv[0]);
        return 1;
    }

    static double noisy_data[MAX_FRAMES][N_MEAS];
    int nframes = 0;
    if (load_csv(argv[1], noisy_data, &nframes) != 0)
        return 1;

    // Small 4x4 blocks
    double F1[16], Q1[16], F1T[16];
    build_F1(F1);
    build_Q1(Q1, PHI);
    asm_mat_transpose(F1, F1T, 4, 4);

    // State and covariance
    double *x = alloc_mat(N_STATE, 1);
    double *P = alloc_mat(N_STATE, N_STATE);
    double *x_pred = alloc_mat(N_STATE, 1);
    double *P_pred = alloc_mat(N_STATE, N_STATE);
    double *x_upd = alloc_mat(N_STATE, 1);
    double *P_upd = alloc_mat(N_STATE, N_STATE);
    build_x0(x, noisy_data[0]);
    build_P0(P);

    // Workspace
    double *Hk = alloc_mat(N_MEAS, N_STATE);
    double *R_adp = alloc_mat(N_MEAS, N_MEAS);
    double *S = alloc_mat(N_MEAS, N_MEAS);
    double *PHt = alloc_mat(N_STATE, N_MEAS);
    double *Kt = alloc_mat(N_MEAS, N_STATE);
    double *z_hat = alloc_mat(N_MEAS, 1);
    double *y = alloc_mat(N_MEAS, 1);
    double *IKH = alloc_mat(N_STATE, N_STATE);
    double *tmp1 = alloc_mat(N_STATE, N_STATE);
    double *tmp2 = alloc_mat(N_STATE, N_STATE);
    double *Ychol = alloc_mat(N_MEAS, N_MEAS);
    double *K = alloc_mat(N_STATE, N_MEAS);
    double z_sph[N_MEAS];
    double temp_3x12[36];

    // Output storage
    double **all_states = (double **)malloc(nframes * sizeof(double *));
    for (int k = 0; k < nframes; k++)
        all_states[k] = alloc_mat(N_STATE, 1);
    mat_copy(all_states[0], x, N_STATE);

    printf("[driver_ekf] Running optimised EKF for %d frames...\n", nframes);
    printf("[driver_ekf] Block-diagonal predict; dense Hk update.\n\n");
    fflush(stdout);

    struct timespec t_start, t_fs, t_now;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    for (int k = 1; k < nframes; k++)
    {
        clock_gettime(CLOCK_MONOTONIC, &t_fs);

        // Predict
        ekf_predict_fast(x, P, F1, Q1, F1T, x_pred, P_pred);

        // Spherical measurement
        cartesian_to_spherical(noisy_data[k], z_sph);

        // Jacobian
        asm_vec_zero_matrix(Hk, (long)N_MEAS * N_STATE);
        asm_jacobian_full(x_pred, Hk, temp_3x12);

        // Adaptive R
        build_R_adaptive(R_adp, x_pred, SIGMA_CART);

        // Update 
        ekf_update(x_pred, P_pred, Hk, R_adp, z_sph,
                   x_upd, P_upd,
                   S, PHt, Kt, z_hat, y, IKH, tmp1, tmp2, Ychol, K);

        mat_copy(x, x_upd, N_STATE);
        mat_copy(P, P_upd, (long)N_STATE * N_STATE);
        mat_copy(all_states[k], x_upd, N_STATE);

        clock_gettime(CLOCK_MONOTONIC, &t_now);
        double fs = (t_now.tv_sec - t_fs.tv_sec) + (t_now.tv_nsec - t_fs.tv_nsec) * 1e-9;
        double el = (t_now.tv_sec - t_start.tv_sec) + (t_now.tv_nsec - t_start.tv_nsec) * 1e-9;
        double eta = fs * (nframes - k);
        printf("  [EKF] frame %4d/%d  frame=%.2fs  elapsed=%5.0fs  ETA=%dm%02ds\n",
               k, nframes - 1, fs, el, (int)(eta / 60), (int)eta % 60);
        fflush(stdout);
    }

    printf("\n[driver_ekf] Done.\n");
    write_output(argv[2], all_states, nframes);

    free(x);
    free(P);
    free(x_pred);
    free(P_pred);
    free(x_upd);
    free(P_upd);
    free(Hk);
    free(R_adp);
    free(S);
    free(PHt);
    free(Kt);
    free(z_hat);
    free(y);
    free(IKH);
    free(tmp1);
    free(tmp2);
    free(Ychol);
    free(K);
    for (int k = 0; k < nframes; k++)
        free(all_states[k]);
    free(all_states);
    return 0;
}