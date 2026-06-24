#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N_JOINTS  23
#define N_AXIS     4
#define N_JOINT   12
#define N_STATE  276
#define N_MEAS    69
#define MAX_FRAMES 4000

#define DT      0.01
#define PHI     0.5
#define R_NOISE 0.05

#define M(A,i,j,cols) ((A)[(long)(i)*(cols)+(j)])

static const char *JOINT_NAMES[N_JOINTS] = {
    "pelvis","L5","L3","T12","T8","neck","head",
    "shoulderRight","upperArmRight","forearmRight","handRight",
    "shoulderLeft","upperArmLeft","forearmLeft","handLeft",
    "upperLegRight","lowerLegRight","footRight","toeRight",
    "upperLegLeft","lowerLegLeft","footLeft","toeLeft"
};

// Assembly function declarations
extern void   asm_mat_mul(double*,double*,double*,long,long,long);
extern void   asm_mat_add(double*,double*,double*,long,long);
extern void   asm_mat_sub(double*,double*,double*,long,long);
extern void   asm_mat_transpose(double*,double*,long,long);
extern void   asm_symmetrize(double*,long);
extern void   asm_cholesky_solve(double*,double*,double*,double*,long,long);
extern double asm_newton_sqrt(double);

// Heap allocator
static double *alloc_mat(long r, long c) {
    double *p = (double*)calloc(r*c, sizeof(double));
    if (!p) { fprintf(stderr,"OOM\n"); exit(1); }
    return p;
}
static void mat_copy(double *dst, const double *src, long n) {
    memcpy(dst, src, n*sizeof(double));
}
static void mat_negate(double *A, long n) {
    for (long i=0;i<n;i++) A[i]=-A[i];
}
static void add_identity(double *A, long n) {
    for (long i=0;i<n;i++) M(A,i,i,n)+=1.0;
}

// Build 4×4 F1 block
static void build_F1(double *F1) {
    double dt=DT,dt2=dt*dt,dt3=dt2*dt;
    memset(F1,0,16*sizeof(double));
    M(F1,0,0,4)=1; M(F1,0,1,4)=dt; M(F1,0,2,4)=dt2/2; M(F1,0,3,4)=dt3/6;
    M(F1,1,1,4)=1; M(F1,1,2,4)=dt; M(F1,1,3,4)=dt2/2;
    M(F1,2,2,4)=1; M(F1,2,3,4)=dt;
    M(F1,3,3,4)=1;
}
static void build_Q1(double *Q1) {
    double dt=DT,dt2=dt*dt,dt3=dt2*dt,dt4=dt3*dt,
           dt5=dt4*dt,dt6=dt5*dt,dt7=dt6*dt,phi=PHI;
    memset(Q1,0,16*sizeof(double));
    M(Q1,0,0,4)=phi*dt7/252; M(Q1,0,1,4)=phi*dt6/72;  M(Q1,1,0,4)=M(Q1,0,1,4);
    M(Q1,0,2,4)=phi*dt5/30;  M(Q1,2,0,4)=M(Q1,0,2,4);
    M(Q1,0,3,4)=phi*dt4/24;  M(Q1,3,0,4)=M(Q1,0,3,4);
    M(Q1,1,1,4)=phi*dt5/20;
    M(Q1,1,2,4)=phi*dt4/8;   M(Q1,2,1,4)=M(Q1,1,2,4);
    M(Q1,1,3,4)=phi*dt3/6;   M(Q1,3,1,4)=M(Q1,1,3,4);
    M(Q1,2,2,4)=phi*dt3/3;   M(Q1,2,3,4)=phi*dt2/2; M(Q1,3,2,4)=M(Q1,2,3,4);
    M(Q1,3,3,4)=phi*dt;
}

// Block-diagonal mat-mul: C = F_blk * B
static void blkdiag_FP(double *F1, double *B, double *C, long cols) {
    // 69 blocks of size 4×4 (one per axis-per-joint)
    for (int b = 0; b < 69; b++) {
        int row = b * 4;       
        double *Bsub = B + row * cols;
        double *Csub = C + row * cols;
        asm_mat_mul(F1, Bsub, Csub, 4, 4, cols);
    }
}

// Block-diagonal P update: P_pred = F*P*F^T + Q
static void predict_P_blockdiag(double *F1, double *Q1,
                                  double *P, double *P_pred,
                                  double *tmp12,   
                                  double *F1T) { 
    
    asm_mat_transpose(F1, F1T, 4, 4);

    // Zero P_pred
    memset(P_pred, 0, (long)N_STATE*N_STATE*sizeof(double));

    double tmp4[16];
    for (int bi = 0; bi < 69; bi++) {
        for (int bj = 0; bj < 69; bj++) {
            int ri = bi*4, rj = bj*4;
            for (int r=0;r<4;r++)
                for (int c=0;c<4;c++)
                    tmp4[r*4+c] = M(P, ri+r, rj+c, N_STATE);

            double fp[16];
            asm_mat_mul(F1, tmp4, fp, 4, 4, 4);

            double fpft[16];
            asm_mat_mul(fp, F1T, fpft, 4, 4, 4);

            for (int r=0;r<4;r++)
                for (int c=0;c<4;c++)
                    M(P_pred, ri+r, rj+c, N_STATE) += fpft[r*4+c];
        }
    }

    for (int b=0;b<69;b++) {
        int base = b*4;
        for (int r=0;r<4;r++)
            for (int c=0;c<4;c++)
                M(P_pred, base+r, base+c, N_STATE) += Q1[r*4+c];
    }

    // Symmetrize
    asm_symmetrize(P_pred, N_STATE);
}

// H*x: exploit block-sparse H 
static void H_times_x(double *x, double *Hx) {
    // H picks px,py,pz from each joint: indices 0,4,8 of each 12-block 
    for (int j=0;j<N_JOINTS;j++) {
        Hx[j*3+0] = x[j*N_JOINT+0];   
        Hx[j*3+1] = x[j*N_JOINT+4];   
        Hx[j*3+2] = x[j*N_JOINT+8];   
    }
}

// P*H^T: result is N×M. H^T col j = sparse 276-vector with 1s at rows
static void P_times_HT(double *P, double *PHt) {
    for (int j=0;j<N_JOINTS;j++) {
        int col_px = j*N_JOINT+0;
        int col_py = j*N_JOINT+4;
        int col_pz = j*N_JOINT+8;
        for (int r=0;r<N_STATE;r++) {
            PHt[r*N_MEAS + j*3+0] = M(P,r,col_px,N_STATE);
            PHt[r*N_MEAS + j*3+1] = M(P,r,col_py,N_STATE);
            PHt[r*N_MEAS + j*3+2] = M(P,r,col_pz,N_STATE);
        }
    }
}

// H*PHt: S = H*PHt + R

static void H_times_PHt(double *PHt, double *R, double *S) {
    for (int i=0;i<N_JOINTS;i++) {
        int src_rows[3] = {i*N_JOINT+0, i*N_JOINT+4, i*N_JOINT+8};
        for (int k=0;k<3;k++) {
            int srow = i*3+k;
            int prow = src_rows[k];
            for (int col=0;col<N_MEAS;col++)
                M(S,srow,col,N_MEAS) = PHt[prow*N_MEAS+col] + M(R,srow,col,N_MEAS);
        }
    }
}

// K*H: N×M × M×N = N×N  (used for IKH = I - K*H)
static void K_times_H(double *K, double *IKH) {

    memset(IKH, 0, (long)N_STATE*N_STATE*sizeof(double));
    for (int j=0;j<N_JOINTS;j++) {
        int meas_cols[3] = {j*3+0, j*3+1, j*3+2};
        int state_cols[3] = {j*N_JOINT+0, j*N_JOINT+4, j*N_JOINT+8};
        for (int r=0;r<N_STATE;r++) {
            for (int k=0;k<3;k++) {
                IKH[r*N_STATE + state_cols[k]] += K[r*N_MEAS + meas_cols[k]];
            }
        }
    }
    // IKH = I - K*H 
    mat_negate(IKH, (long)N_STATE*N_STATE);
    add_identity(IKH, N_STATE);
}

// Build R (diagonal, constant) 
static void build_R(double *R) {
    memset(R, 0, (long)N_MEAS*N_MEAS*sizeof(double));
    for (int i=0;i<N_MEAS;i++) M(R,i,i,N_MEAS) = R_NOISE;
}

static void build_P0(double *P) {
    memset(P, 0, (long)N_STATE*N_STATE*sizeof(double));
    for (int j=0;j<N_JOINTS;j++) {
        int b=j*N_JOINT;
        for (int a=0;a<3;a++) {
            int ab=b+a*N_AXIS;
            M(P,ab+0,ab+0,N_STATE)=1.0;
            M(P,ab+1,ab+1,N_STATE)=100.0;
            M(P,ab+2,ab+2,N_STATE)=100.0;
            M(P,ab+3,ab+3,N_STATE)=100.0;
        }
    }
}
static void build_x0(double *x, double *z0) {
    memset(x,0,N_STATE*sizeof(double));
    for (int j=0;j<N_JOINTS;j++) {
        x[j*N_JOINT+0]=z0[j*3+0];
        x[j*N_JOINT+4]=z0[j*3+1];
        x[j*N_JOINT+8]=z0[j*3+2];
    }
}

// Fast LKF predict
static void lkf_predict_fast(double *x, double *P,
                               double *F1, double *Q1,
                               double *x_pred, double *P_pred,
                               double *F1T) {
    for (int b=0;b<69;b++) {
        int row=b*4;
        double xblk[4]={x[row],x[row+1],x[row+2],x[row+3]};
        double yblk[4];
        asm_mat_mul(F1, xblk, yblk, 4, 4, 1);
        x_pred[row]=yblk[0]; x_pred[row+1]=yblk[1];
        x_pred[row+2]=yblk[2]; x_pred[row+3]=yblk[3];
    }

    // P_pred = F*P*F^T + Q (block structure)
    predict_P_blockdiag(F1, Q1, P, P_pred, NULL, F1T);
}

// Fast LKF update
static void lkf_update_fast(double *x_pred, double *P_pred,
                              double *R, double *z,
                              double *x_upd, double *P_upd,
                              double *PHt, double *S, double *Ychol,
                              double *Kt, double *K,
                              double *Hxp, double *IKH,
                              double *tmp1, double *tmp2) {
    // PHt = P * H^T  (exploit H sparsity)
    P_times_HT(P_pred, PHt);

    // S = H * PHt + R  (exploit H sparsity)
    H_times_PHt(PHt, R, S);
    asm_symmetrize(S, N_MEAS);

    // Solve S * K^T = PHt^T 
    asm_mat_transpose(PHt, Kt, N_STATE, N_MEAS); 
    mat_copy(Ychol, S, (long)N_MEAS*N_MEAS);
    asm_cholesky_solve(Ychol, Kt, Kt, tmp1, N_MEAS, N_STATE);
    asm_mat_transpose(Kt, K, N_MEAS, N_STATE);  

    // y = z - H*x_pred  (exploit H sparsity)
    H_times_x(x_pred, Hxp);
    asm_mat_sub(z, Hxp, Hxp, N_MEAS, 1);

    // x_upd = x_pred + K*y  
    asm_mat_mul(K, Hxp, x_upd, N_STATE, N_MEAS, 1);
    asm_mat_add(x_pred, x_upd, x_upd, N_STATE, 1);

    // IKH = I - K*H  (exploit H sparsity)
    K_times_H(K, IKH);

    // P_upd = IKH * P_pred * IKH^T + K*R*K^T 
    asm_mat_mul(IKH, P_pred, tmp1, N_STATE, N_STATE, N_STATE);
    asm_mat_transpose(IKH, tmp2, N_STATE, N_STATE);
    asm_mat_mul(tmp1, tmp2, P_upd, N_STATE, N_STATE, N_STATE);

    // K*R*K^T: R is diagonal so K*R = K scaled per column
    memset(tmp1, 0, (long)N_STATE*N_STATE*sizeof(double));
    for (int i=0;i<N_MEAS;i++) {
        double rii = M(R,i,i,N_MEAS);
        for (int r=0;r<N_STATE;r++)
            for (int c=0;c<N_STATE;c++)
                tmp1[r*N_STATE+c] += rii * K[r*N_MEAS+i] * K[c*N_MEAS+i];
    }
    asm_mat_add(P_upd, tmp1, P_upd, N_STATE, N_STATE);
    asm_symmetrize(P_upd, N_STATE);
}

// CSV load 
static int load_csv(const char *fn, double data[][N_MEAS], int *nf) {
    FILE *f=fopen(fn,"r"); if(!f){perror(fn);return -1;}
    char line[8192]; int frame=0,hdr=0;
    while(fgets(line,sizeof(line),f)){
        if(!hdr){char*p=line;while(*p==' '||*p=='\t')p++;
            if(*p<'0'||*p>'9'){hdr=1;continue;}hdr=1;}
        if(frame>=MAX_FRAMES)break;
        char*tok=strtok(line,",\n\r");
        for(int i=0;i<N_MEAS&&tok;i++){data[frame][i]=atof(tok);tok=strtok(NULL,",\n\r");}
        frame++;
    }
    fclose(f);*nf=frame;
    printf("[driver_lkf] Loaded %d frames from '%s'\n",frame,fn);
    return 0;
}

// CSV write
static void write_output(const char *fn, double **states, int nf) {
    FILE *f=fopen(fn,"w");if(!f){perror(fn);return;}
    fprintf(f,"frame");
    const char *ax[]={"x","y","z"},*st[]={"p","v","a","j"};
    for(int j=0;j<N_JOINTS;j++)
        for(int a=0;a<3;a++)
            for(int s=0;s<N_AXIS;s++)
                fprintf(f,",%s_%s%s",JOINT_NAMES[j],st[s],ax[a]);
    fprintf(f,"\n");
    for(int k=0;k<nf;k++){
        fprintf(f,"%d",k);
        for(int i=0;i<N_STATE;i++) fprintf(f,",%.8f",states[k][i]);
        fprintf(f,"\n");
    }
    fclose(f);
    printf("[driver_lkf] Output written to '%s'\n",fn);
}


int main(int argc, char *argv[]) {
    printf("=========================================\n");
    printf("  Team Thunder — LKF Assembly Driver\n");
    printf("  Milestone 3 — Block-Diagonal Optimised\n");
    printf("=========================================\n\n");
    if(argc<3){fprintf(stderr,"Usage: %s <noisy.csv> <out.csv>\n",argv[0]);return 1;}

    static double noisy_data[MAX_FRAMES][N_MEAS];
    int nframes=0;
    if(load_csv(argv[1],noisy_data,&nframes)!=0) return 1;

    // Build small (4×4) F1 and Q1 blocks
    double F1[16], Q1[16], F1T[16];
    build_F1(F1);
    build_Q1(Q1);
    asm_mat_transpose(F1, F1T, 4, 4);

    // Build R (69×69 diagonal)
    double *R = alloc_mat(N_MEAS, N_MEAS);
    build_R(R);

    // State and covariance
    double *x     = alloc_mat(N_STATE, 1);
    double *P     = alloc_mat(N_STATE, N_STATE);
    double *x_pred= alloc_mat(N_STATE, 1);
    double *P_pred= alloc_mat(N_STATE, N_STATE);
    double *x_upd = alloc_mat(N_STATE, 1);
    double *P_upd = alloc_mat(N_STATE, N_STATE);
    build_x0(x, noisy_data[0]);
    build_P0(P);

    // Workspace
    double *PHt   = alloc_mat(N_STATE, N_MEAS);
    double *S     = alloc_mat(N_MEAS,  N_MEAS);
    double *Ychol = alloc_mat(N_MEAS,  N_MEAS);
    double *Kt    = alloc_mat(N_MEAS,  N_STATE);
    double *K     = alloc_mat(N_STATE, N_MEAS);
    double *Hxp   = alloc_mat(N_MEAS,  1);
    double *IKH   = alloc_mat(N_STATE, N_STATE);
    double *tmp1  = alloc_mat(N_STATE, N_STATE);
    double *tmp2  = alloc_mat(N_STATE, N_STATE);

    // Output storage
    double **all_states=(double**)malloc(nframes*sizeof(double*));
    for(int k=0;k<nframes;k++) all_states[k]=alloc_mat(N_STATE,1);
    mat_copy(all_states[0], x, N_STATE);

    printf("[driver_lkf] Running optimised LKF for %d frames...\n",nframes);
    printf("[driver_lkf] Block-diagonal predict + sparse H update.\n\n");
    fflush(stdout);

    struct timespec t_start, t_fs, t_now;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    for(int k=1;k<nframes;k++){
        clock_gettime(CLOCK_MONOTONIC, &t_fs);

        lkf_predict_fast(x, P, F1, Q1, x_pred, P_pred, F1T);

        double *z = noisy_data[k];

        lkf_update_fast(x_pred, P_pred, R, z,
                        x_upd, P_upd,
                        PHt, S, Ychol, Kt, K, Hxp, IKH, tmp1, tmp2);

        mat_copy(x,     x_upd, N_STATE);
        mat_copy(P,     P_upd, (long)N_STATE*N_STATE);
        mat_copy(all_states[k], x_upd, N_STATE);

        clock_gettime(CLOCK_MONOTONIC, &t_now);
        double fs  = (t_now.tv_sec - t_fs.tv_sec) + (t_now.tv_nsec - t_fs.tv_nsec) * 1e-9;
        double el  = (t_now.tv_sec - t_start.tv_sec) + (t_now.tv_nsec - t_start.tv_nsec) * 1e-9;
        double eta = fs * (nframes - k - 1);
        
        int elapsed_min = (int)(el / 60);
        int elapsed_sec = (int)el % 60;
        int eta_min = (int)(eta / 60);
        int eta_sec = (int)eta % 60;
        
        printf("  [LKF] frame %4d/%d  frame=%.2fs  elapsed=%02d:%02d  ETA=%02d:%02d\n",
               k, nframes-1, fs, elapsed_min, elapsed_sec, eta_min, eta_sec);
        fflush(stdout);
    }

    printf("\n[driver_lkf] Done.\n");
    write_output(argv[2], all_states, nframes);

    free(R);free(x);free(P);free(x_pred);free(P_pred);free(x_upd);free(P_upd);
    free(PHt);free(S);free(Ychol);free(Kt);free(K);free(Hxp);
    free(IKH);free(tmp1);free(tmp2);
    for(int k=0;k<nframes;k++) free(all_states[k]);
    free(all_states);
    return 0;
}