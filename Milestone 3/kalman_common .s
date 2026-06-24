    .section .text
    .global asm_mat_mul
asm_mat_mul:
    addi    sp, sp, -96
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1,  16(sp)
    sd      s2,  24(sp)
    sd      s3,  32(sp)
    sd      s4,  40(sp)
    sd      s5,  48(sp)
    sd      s6,  56(sp)
    sd      s7,  64(sp)
    sd      s8,  72(sp)
    sd      s9,  80(sp)
    sd      s10, 88(sp)

    mv      s0, a0          
    mv      s1, a1          
    mv      s2, a2         
    mv      s3, a3          
    mv      s4, a4          
    mv      s5, a5          

    # Zero C (M*N doubles) using pointer walk
    mul     t0, s3, s5      
    slli    t0, t0, 3       
    add     t0, s2, t0     
    mv      t1, s2
    fmv.d.x ft3, zero
.zero_C:
    bge     t1, t0, .zero_done
    fsd     ft3, 0(t1)
    addi    t1, t1, 8
    j       .zero_C
.zero_done:

    # i loop: rows of A and C
    li      s6, 0           
    mv      s9,  s0         
    mv      s10, s2         
    slli    t6, s5, 3     

.loop_i:
    bge     s6, s3, .done_mat_mul

    # k loop: shared dimension 
    li      s8, 0          
    mv      t3, s1          

.loop_k:
    bge     s8, s4, .next_i

    # Load A[i][k]: A_row + k*8
    slli    t0, s8, 3
    add     t0, s9, t0
    fld     ft1, 0(t0)      

    # j loop: cols of B (pointer walk, zero multiplies)
    mv      t4, t3          
    mv      t5, s10        
    li      s7, 0

.loop_j:
    bge     s7, s5, .next_k

    fld     ft2, 0(t4)      
    fld     ft0, 0(t5)      
    fmadd.d ft0, ft1, ft2, ft0   
    fsd     ft0, 0(t5)

    addi    t4, t4, 8       
    addi    t5, t5, 8       
    addi    s7, s7, 1
    j       .loop_j

.next_k:
    # t6 holds N*8
    add     t3, t3, t6      
    addi    s8, s8, 1
    j       .loop_k

.next_i:
    slli    t0, s4, 3       
    add     s9,  s9,  t0    
    add     s10, s10, t6    
    addi    s6, s6, 1
    j       .loop_i

.done_mat_mul:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
    ld      s1,  16(sp)
    ld      s2,  24(sp)
    ld      s3,  32(sp)
    ld      s4,  40(sp)
    ld      s5,  48(sp)
    ld      s6,  56(sp)
    ld      s7,  64(sp)
    ld      s8,  72(sp)
    ld      s9,  80(sp)
    ld      s10, 88(sp)
    addi    sp, sp, 96
    ret

    .global asm_mat_add
asm_mat_add:
    # Pointer-walk: no index multiply, just advance 3 pointers by 8 each step
    mul     t0, a3, a4      
    slli    t0, t0, 3       
    add     t0, a2, t0     
    mv      t1, a0          
    mv      t2, a1          
    mv      t3, a2          
.loop_add:
    bge     t3, t0, .done_add
    fld     ft0, 0(t1)
    fld     ft1, 0(t2)
    fadd.d  ft0, ft0, ft1
    fsd     ft0, 0(t3)
    addi    t1, t1, 8
    addi    t2, t2, 8
    addi    t3, t3, 8
    j       .loop_add
.done_add:
    ret

    .global asm_mat_sub
asm_mat_sub:
    mul     t0, a3, a4
    slli    t0, t0, 3
    add     t0, a2, t0      # end pointer
    mv      t1, a0
    mv      t2, a1
    mv      t3, a2
.loop_sub:
    bge     t3, t0, .done_sub
    fld     ft0, 0(t1)
    fld     ft1, 0(t2)
    fsub.d  ft0, ft0, ft1
    fsd     ft0, 0(t3)
    addi    t1, t1, 8
    addi    t2, t2, 8
    addi    t3, t3, 8
    j       .loop_sub
.done_sub:
    ret

    .global asm_mat_transpose
asm_mat_transpose:
    # Optimised: A_ptr walks row i left-to-right (stride 8).
    # B write address = B + (j*M + i)*8; j increments so B stride = M*8.
    # We recompute B write base once per i (B + i*8), then stride by M*8 per j.
    addi    sp, sp, -48
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)

    mv      s0, a0          
    mv      s1, a1         
    mv      s2, a2          
    mv      s3, a3          
    slli    t6, s2, 3       

    li      s4, 0          
.loop_tp_i:
    bge     s4, s2, .done_tp

    # A_ptr = A + i*N*8  (start of row i)
    mul     t0, s4, s3
    slli    t0, t0, 3
    add     t1, s0, t0      

    # B_ptr = B + i*8  (column i of B, first element B[0][i])
    slli    t0, s4, 3
    add     t2, s1, t0      

    li      t3, 0          
.loop_tp_j:
    bge     t3, s3, .next_tp_i
    fld     ft0, 0(t1)     
    fsd     ft0, 0(t2)     
    addi    t1, t1, 8      
    add     t2, t2, t6    
    addi    t3, t3, 1
    j       .loop_tp_j
.next_tp_i:
    addi    s4, s4, 1
    j       .loop_tp_i
.done_tp:
    ld      ra, 0(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    addi    sp, sp, 48
    ret

    .global asm_symmetrize
asm_symmetrize:
    addi    sp, sp, -32
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)

    mv      s0, a0
    mv      s1, a1

    # Load 0.5 constant
    li      t0, 0
    la      t1, .const_half
    fld     fs0, 0(t1)      

    li      s2, 0           
.sym_i:
    bge     s2, s1, .sym_done
    addi    t0, s2, 1       
.sym_j:
    bge     t0, s1, .sym_next_i
    # A[i][j] offset
    mul     t1, s2, s1
    add     t1, t1, t0
    slli    t1, t1, 3
    add     t1, s0, t1
    # A[j][i] offset
    mul     t2, t0, s1
    add     t2, t2, s2
    slli    t2, t2, 3
    add     t2, s0, t2

    fld     ft0, 0(t1)      
    fld     ft1, 0(t2)      
    fadd.d  ft0, ft0, ft1   
    fmul.d  ft0, ft0, fs0  
    fsd     ft0, 0(t1)
    fsd     ft0, 0(t2)

    addi    t0, t0, 1
    j       .sym_j
.sym_next_i:
    addi    s2, s2, 1
    j       .sym_i
.sym_done:
    ld      ra, 0(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    addi    sp, sp, 32
    ret

    .global asm_newton_sqrt
asm_newton_sqrt:
    # Load constants
    la      t0, .const_half
    fld     fa3, 0(t0)          

    # Check for zero input
    la      t0, .const_zero
    fld     ft0, 0(t0)
    feq.d   t1, fa0, ft0        
    bnez    t1, .sqrt_zero

    # Initial guess: x = s / 2
    fmul.d  fa1, fa0, fa3       

    li      t0, 0               
    li      t1, 20              

.sqrt_iter:
    bge     t0, t1, .sqrt_done
    # x_new = 0.5 * (x + s/x)
    fdiv.d  fa4, fa0, fa1       
    fadd.d  fa2, fa1, fa4       
    fmul.d  fa2, fa2, fa3       

    # Check convergence: |x_new - x| < 1e-15 * x
    fsub.d  ft1, fa2, fa1      
    # Take absolute value of ft1
    la      t2, .const_neg_one
    fld     ft2, 0(t2)
    flt.d   t3, ft1, ft0        
    beqz    t3, .sqrt_abs_done
    fmul.d  ft1, ft1, ft2       
.sqrt_abs_done:
    # tolerance = 1e-15 * fa1 (relative)
    la      t2, .const_eps
    fld     ft2, 0(t2)
    fmul.d  ft2, ft2, fa1      
    flt.d   t3, ft1, ft2        
    fmv.d   fa1, fa2            
    bnez    t3, .sqrt_done      
    addi    t0, t0, 1
    j       .sqrt_iter

.sqrt_done:
    # fa0 still holds s; return value in fa0
    fmv.d   fa0, fa1           
    ret

.sqrt_zero:
    # Return 0.0
    la      t0, .const_zero
    fld     fa0, 0(t0)
    ret

    .global asm_cholesky
asm_cholesky:
    addi    sp, sp, -64
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)
    sd      s5, 48(sp)

    mv      s0, a0          
    mv      s1, a1         

    la      t0, .const_half
    fld     fs0, 0(t0)      

    li      s2, 0       
.chol_j:
    bge     s2, s1, .chol_done

    # Compute diagonal element L[j][j] 
    # A[j][j] address
    mul     t0, s2, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fld     ft0, 0(t0)      

    # Subtract sum_{k=0}^{j-1} L[j][k]^2
    li      s3, 0           # k = 0
.chol_diag_k:
    bge     s3, s2, .chol_diag_done
    mul     t1, s2, s1
    add     t1, t1, s3
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)      
    # fmsub.d: ft0 = ft1 * ft1 - ft0 
    # Use fmadd.d with negation: ft0 = ft1*ft1 + ft0 would add, so:
    fmul.d  ft2, ft1, ft1   
    fsub.d  ft0, ft0, ft2  
    addi    s3, s3, 1
    j       .chol_diag_k
.chol_diag_done:
    # ft0 = A[j][j] - sum  => must be > 0 for SPD
    # L[j][j] = sqrt(ft0)
    fmv.d   fa0, ft0
    call    asm_newton_sqrt   
    fmv.d   ft3, fa0        

    # Store L[j][j]
    mul     t0, s2, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fsd     ft3, 0(t0)

    # Compute off-diagonal elements L[i][j] for i > j 
    addi    s4, s2, 1       # i = j + 1
.chol_i:
    bge     s4, s1, .chol_next_j

    # Start with A[i][j]
    mul     t0, s4, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fld     ft4, 0(t0)    

    # Subtract sum_{k=0}^{j-1} L[i][k] * L[j][k]
    li      s5, 0           
.chol_off_k:
    bge     s5, s2, .chol_off_done
    # L[i][k]
    mul     t1, s4, s1
    add     t1, t1, s5
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft5, 0(t1)
    # L[j][k]
    mul     t2, s2, s1
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s0, t2
    fld     ft6, 0(t2)
    # ft4 -= L[i][k] * L[j][k]
    # Use fnmsub.d: fd = -(fs1*fs2) + fs3  which gives fs3 - fs1*fs2
    fnmsub.d ft4, ft5, ft6, ft4   
    addi    s5, s5, 1
    j       .chol_off_k
.chol_off_done:
    # L[i][j] = ft4 / L[j][j]
    fdiv.d  ft4, ft4, ft3
    # Store L[i][j]
    mul     t0, s4, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fsd     ft4, 0(t0)

    addi    s4, s4, 1
    j       .chol_i

.chol_next_j:
    addi    s2, s2, 1
    j       .chol_j

.chol_done:
    ld      ra, 0(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    ld      s5, 48(sp)
    addi    sp, sp, 64
    ret

    .global asm_forward_sub
asm_forward_sub:
    addi    sp, sp, -56
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)
    sd      s5, 48(sp)

    mv      s0, a0         
    mv      s1, a1         
    mv      s2, a2          
    mv      s3, a3        
    mv      s4, a4         

    li      s5, 0           
.fsub_col:
    bge     s5, s4, .fsub_done

    li      t6, 0          
.fsub_row:
    bge     t6, s3, .fsub_next_col

    # Load B[i][c]
    mul     t0, t6, s4
    add     t0, t0, s5
    slli    t0, t0, 3
    add     t0, s1, t0
    fld     ft0, 0(t0)      

    # Subtract sum_{k=0}^{i-1} L[i][k] * Y[k][c]
    li      t5, 0          
.fsub_k:
    bge     t5, t6, .fsub_k_done
    # L[i][k]
    mul     t1, t6, s3
    add     t1, t1, t5
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    # Y[k][c]
    mul     t2, t5, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fld     ft2, 0(t2)
    # ft0 -= L[i][k] * Y[k][c]
    fnmsub.d ft0, ft1, ft2, ft0
    addi    t5, t5, 1
    j       .fsub_k
.fsub_k_done:
    # Divide by L[i][i]
    mul     t1, t6, s3
    add     t1, t1, t6
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)      
    fdiv.d  ft0, ft0, ft1
    # Store Y[i][c]
    mul     t2, t6, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fsd     ft0, 0(t2)

    addi    t6, t6, 1
    j       .fsub_row
.fsub_next_col:
    addi    s5, s5, 1
    j       .fsub_col
.fsub_done:
    ld      ra, 0(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    ld      s5, 48(sp)
    addi    sp, sp, 56
    ret

    .global asm_back_sub
asm_back_sub:
    addi    sp, sp, -56
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)
    sd      s5, 48(sp)

    mv      s0, a0          
    mv      s1, a1          
    mv      s2, a2          
    mv      s3, a3          
    mv      s4, a4         

    li      s5, 0           
.bsub_col:
    bge     s5, s4, .bsub_done

    addi    t6, s3, -1      
.bsub_row:
    li      t0, 0
    blt     t6, t0, .bsub_next_col  

    # Load Y[i][c]
    mul     t0, t6, s4
    add     t0, t0, s5
    slli    t0, t0, 3
    add     t0, s1, t0
    fld     ft0, 0(t0)

    # Subtract sum_{k=i+1}^{N-1} L[k][i] * X[k][c]
    addi    t5, t6, 1      
.bsub_k:
    bge     t5, s3, .bsub_k_done
    # L[k][i] = L[k*N + i]  (i-th col of L row k)
    mul     t1, t5, s3
    add     t1, t1, t6
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    # X[k][c]
    mul     t2, t5, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fld     ft2, 0(t2)
    fnmsub.d ft0, ft1, ft2, ft0
    addi    t5, t5, 1
    j       .bsub_k
.bsub_k_done:
    # Divide by L[i][i]
    mul     t1, t6, s3
    add     t1, t1, t6
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    fdiv.d  ft0, ft0, ft1
    # Store X[i][c]
    mul     t2, t6, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fsd     ft0, 0(t2)

    addi    t6, t6, -1      
    j       .bsub_row
.bsub_next_col:
    addi    s5, s5, 1
    j       .bsub_col
.bsub_done:
    ld      ra, 0(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    ld      s5, 48(sp)
    addi    sp, sp, 56
    ret

    .global asm_cholesky_solve
asm_cholesky_solve:
    addi    sp, sp, -48
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)

    mv      s0, a0          
    mv      s1, a1          
    mv      s2, a2         
    mv      s3, a3          
    mv      s4, a4          
    mv      s5, a5          

    # Step 1: Cholesky decomposition of A 
    mv      a0, s0
    mv      a1, s4
    call    asm_cholesky

    # Step 2: Forward substitution — solve L * Y = B
    mv      a0, s0          
    mv      a1, s1          
    mv      a2, s3          
    mv      a3, s4          
    mv      a4, s5          
    call    asm_forward_sub

    # Step 3: Backward substitution — solve L^T * X = Y
    mv      a0, s0         
    mv      a1, s3          
    mv      a2, s2          
    mv      a3, s4        
    mv      a4, s5          
    call    asm_back_sub

    ld      ra, 0(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    addi    sp, sp, 48
    ret

# Read-only constants used by multiple functions
    .section .rodata
    .align 3                  

.const_half:
    .double 0.5

.const_zero:
    .double 0.0

.const_neg_one:
    .double -1.0

.const_eps:
    .double 1e-15

.const_two:
    .double 2.0