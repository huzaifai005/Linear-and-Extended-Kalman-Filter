    .section .text


    .global asm_mat_mul
asm_mat_mul:
    addi    sp, sp, -96
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)
    sd      s5, 48(sp)
    sd      s6, 56(sp)
    sd      s7, 64(sp)
    sd      s8, 72(sp)
    sd      s9, 80(sp)
    sd      s10,88(sp)

    mv      s0, a0          # A
    mv      s1, a1          # B
    mv      s2, a2          # C
    mv      s3, a3          # M
    mv      s4, a4          # K
    mv      s5, a5          # N

    slli    s10, s5, 3      # N*8 (byte stride for B rows and C rows)
    slli    s11, s4, 3      # K*8 (byte stride for A rows)

    # Zero C using vector stores
    mul     t0, s3, s5      # total elements = M * N
    mv      t1, s2          # pC = C
    fmv.d.x ft3, zero       # ft3 = 0.0
.vec_zero_C:
    beqz    t0, .vec_zero_done
    vsetvli t2, t0, e64, m4, ta, ma  
    vfmv.v.f v0, ft3           
    vse64.v v0, (t1)
    slli    t3, t2, 3
    add     t1, t1, t3      # advance pC
    sub     t0, t0, t2      # remaining elements
    j       .vec_zero_C
.vec_zero_done:

    # Outer i loop
    li      s6, 0           # i = 0
    mv      s8, s0          # A_row = A
    mv      s9, s2          # C_row = C

.vec_loop_i:
    bge     s6, s3, .vec_done_matmul

    # Middle k loop 
    li      s7, 0           # k = 0
    mv      t0, s1          # B_row = B

.vec_loop_k:
    bge     s7, s4, .vec_next_i

    # Load scalar A[i][k]
    slli    t3, s7, 3       # k*8
    add     t3, s8, t3
    fld     ft0, 0(t3)      # ft0 = A[i][k]

    # Inner j loop: vectorised
    # Process entire row k of B, accumulating into row i of C.
    # Each iteration handles VL doubles.
    mv      t1, s5        
    mv      t3, s9          
    mv      t4, t0         

.vec_loop_j:
    beqz    t1, .vec_next_k
    vsetvli t2, t1, e64, m4, ta, ma   # m4: 4 regs merged = 4x VL

    vle64.v v0, (t4)        
    vle64.v v4, (t3)        
    vfmacc.vf v4, ft0, v0  
    vse64.v v4, (t3)        

    slli    t5, t2, 3       
    add     t3, t3, t5      
    add     t4, t4, t5      
    sub     t1, t1, t2      
    j       .vec_loop_j

.vec_next_k:
    add     t0, t0, s10     # B_row += N*8
    addi    s7, s7, 1
    j       .vec_loop_k

.vec_next_i:
    add     s8, s8, s11     # A_row += K*8
    add     s9, s9, s10     # C_row += N*8
    addi    s6, s6, 1
    j       .vec_loop_i

.vec_done_matmul:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    ld      s5, 48(sp)
    ld      s6, 56(sp)
    ld      s7, 64(sp)
    ld      s8, 72(sp)
    ld      s9, 80(sp)
    ld      s10,88(sp)
    addi    sp, sp, 96
    ret


    .global asm_mat_add
asm_mat_add:
    mul     t0, a3, a4      # total elements
    mv      t1, a0          # pA
    mv      t2, a1          # pB
    mv      t3, a2          # pC
.vec_loop_add:
    beqz    t0, .vec_done_add
    vsetvli t4, t0, e64, m4, ta, ma
    vle64.v v0, (t1)
    vle64.v v4, (t2)
    vfadd.vv v8, v0, v4
    vse64.v v8, (t3)
    slli    t5, t4, 3
    add     t1, t1, t5
    add     t2, t2, t5
    add     t3, t3, t5
    sub     t0, t0, t4
    j       .vec_loop_add
.vec_done_add:
    ret


    .global asm_mat_sub
asm_mat_sub:
    mul     t0, a3, a4
    mv      t1, a0
    mv      t2, a1
    mv      t3, a2
.vec_loop_sub:
    beqz    t0, .vec_done_sub
    vsetvli t4, t0, e64, m4, ta, ma
    vle64.v v0, (t1)
    vle64.v v4, (t2)
    vfsub.vv v8, v0, v4
    vse64.v v8, (t3)
    slli    t5, t4, 3
    add     t1, t1, t5
    add     t2, t2, t5
    add     t3, t3, t5
    sub     t0, t0, t4
    j       .vec_loop_sub
.vec_done_sub:
    ret


    .global asm_mat_transpose
asm_mat_transpose:
    addi    sp, sp, -48
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)

    mv      s0, a0          # A
    mv      s1, a1          # B
    mv      s2, a2          # M (rows of A)
    mv      s3, a3          # N (cols of A)
    slli    s5, s2, 3      

    li      s4, 0           

.vec_tp_i:
    bge     s4, s2, .vec_tp_done

    # A_ptr = A + i*N*8
    mul     t0, s4, s3
    slli    t0, t0, 3
    add     t0, s0, t0      

    # B_ptr = B + i*8  (start of column i in B)
    slli    t1, s4, 3
    add     t1, s1, t1    

    mv      t2, s3        

.vec_tp_j:
    beqz    t2, .vec_tp_next_i
    vsetvli t3, t2, e64, m1, ta, ma    # t3 = vl

    vle64.v  v0, (t0)       # load vl elements of A row i (unit stride)
    vsse64.v v0, (t1), s5   # scatter into column i of B (stride = M*8)

    slli    t4, t3, 3       
    add     t0, t0, t4      

    mul     t5, t3, s5    
    add     t1, t1, t5    

    sub     t2, t2, t3
    j       .vec_tp_j

.vec_tp_next_i:
    addi    s4, s4, 1
    j       .vec_tp_i

.vec_tp_done:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    addi    sp, sp, 48
    ret


    .global asm_symmetrize
asm_symmetrize:
    addi    sp, sp, -48
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)

    mv      s0, a0
    mv      s1, a1
    slli    s3, s1, 3       # s3 = N*8 (column stride)

    # Load 0.5 into fs0
    la      t0, .vec_half
    fld     fs0, 0(t0)

    li      s2, 0           # i = 0

.vec_sym_i:
    bge     s2, s1, .vec_sym_done

    # Row pointer: A[i][i+1] = A + i*N*8 + (i+1)*8
    mul     t0, s2, s1     
    add     t0, t0, s2      
    addi    t0, t0, 1      
    slli    t0, t0, 3
    add     t0, s0, t0     

    # Column pointer: A[i+1][i] = A + (i+1)*N*8 + i*8
    addi    t6, s2, 1     
    mul     t1, t6, s1     
    add     t1, t1, s2      
    slli    t1, t1, 3
    add     t1, s0, t1      

    # remaining elements = N - i - 1
    sub     t2, s1, s2
    addi    t2, t2, -1

.vec_sym_j:
    beqz    t2, .vec_sym_next_i
    vsetvli t3, t2, e64, m1, ta, ma

    vle64.v  v0, (t0)       
    vlse64.v v1, (t1), s3   

    vfadd.vv v2, v0, v1     
    vfmul.vf v2, v2, fs0    

    vse64.v  v2, (t0)      
    vsse64.v v2, (t1), s3   

    slli    t4, t3, 3       
    add     t0, t0, t4     

    mul     t5, t3, s3     
    add     t1, t1, t5      

    sub     t2, t2, t3
    j       .vec_sym_j

.vec_sym_next_i:
    addi    s2, s2, 1
    j       .vec_sym_i

.vec_sym_done:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    addi    sp, sp, 48
    ret


    .global asm_newton_sqrt
asm_newton_sqrt:
    la      t0, .vec_half
    fld     fa3, 0(t0)
    la      t0, .vec_zero
    fld     ft0, 0(t0)
    feq.d   t1, fa0, ft0
    bnez    t1, .vsqrt_zero
    fmul.d  fa1, fa0, fa3
    li      t0, 0
    li      t1, 20
.vsqrt_iter:
    bge     t0, t1, .vsqrt_done
    fdiv.d  fa4, fa0, fa1
    fadd.d  fa2, fa1, fa4
    fmul.d  fa2, fa2, fa3
    fsub.d  ft1, fa2, fa1
    la      t2, .vec_neg_one
    fld     ft2, 0(t2)
    flt.d   t3, ft1, ft0
    beqz    t3, .vsqrt_abs_done
    fmul.d  ft1, ft1, ft2
.vsqrt_abs_done:
    la      t2, .vec_eps
    fld     ft2, 0(t2)
    fmul.d  ft2, ft2, fa1
    flt.d   t3, ft1, ft2
    fmv.d   fa1, fa2
    bnez    t3, .vsqrt_done
    addi    t0, t0, 1
    j       .vsqrt_iter
.vsqrt_done:
    fmv.d   fa0, fa1
    ret
.vsqrt_zero:
    la      t0, .vec_zero
    fld     fa0, 0(t0)
    ret

    .global asm_cholesky
asm_cholesky:
    addi    sp, sp, -64
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)
    sd      s5, 48(sp)
    mv      s0, a0
    mv      s1, a1
    li      s2, 0
.vchol_j:
    bge     s2, s1, .vchol_done
    mul     t0, s2, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fld     ft0, 0(t0)
    li      s3, 0
.vchol_diag_k:
    bge     s3, s2, .vchol_diag_done
    mul     t1, s2, s1
    add     t1, t1, s3
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    fmul.d  ft2, ft1, ft1
    fsub.d  ft0, ft0, ft2
    addi    s3, s3, 1
    j       .vchol_diag_k
.vchol_diag_done:
    fmv.d   fa0, ft0
    call    asm_newton_sqrt
    fmv.d   ft3, fa0
    mul     t0, s2, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fsd     ft3, 0(t0)
    addi    s4, s2, 1
.vchol_i:
    bge     s4, s1, .vchol_next_j
    mul     t0, s4, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fld     ft4, 0(t0)
    li      s5, 0
.vchol_off_k:
    bge     s5, s2, .vchol_off_done
    mul     t1, s4, s1
    add     t1, t1, s5
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft5, 0(t1)
    mul     t2, s2, s1
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s0, t2
    fld     ft6, 0(t2)
    fnmsub.d ft4, ft5, ft6, ft4
    addi    s5, s5, 1
    j       .vchol_off_k
.vchol_off_done:
    fdiv.d  ft4, ft4, ft3
    mul     t0, s4, s1
    add     t0, t0, s2
    slli    t0, t0, 3
    add     t0, s0, t0
    fsd     ft4, 0(t0)
    addi    s4, s4, 1
    j       .vchol_i
.vchol_next_j:
    addi    s2, s2, 1
    j       .vchol_j
.vchol_done:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
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
    sd      ra,  0(sp)
    sd      s0,  8(sp)
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
.vfsub_col:
    bge     s5, s4, .vfsub_done
    li      t6, 0
.vfsub_row:
    bge     t6, s3, .vfsub_next_col
    mul     t0, t6, s4
    add     t0, t0, s5
    slli    t0, t0, 3
    add     t0, s1, t0
    fld     ft0, 0(t0)
    li      t5, 0
.vfsub_k:
    bge     t5, t6, .vfsub_k_done
    mul     t1, t6, s3
    add     t1, t1, t5
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    mul     t2, t5, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fld     ft2, 0(t2)
    fnmsub.d ft0, ft1, ft2, ft0
    addi    t5, t5, 1
    j       .vfsub_k
.vfsub_k_done:
    mul     t1, t6, s3
    add     t1, t1, t6
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    fdiv.d  ft0, ft0, ft1
    mul     t2, t6, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fsd     ft0, 0(t2)
    addi    t6, t6, 1
    j       .vfsub_row
.vfsub_next_col:
    addi    s5, s5, 1
    j       .vfsub_col
.vfsub_done:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
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
    sd      ra,  0(sp)
    sd      s0,  8(sp)
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
.vbsub_col:
    bge     s5, s4, .vbsub_done
    addi    t6, s3, -1
.vbsub_row:
    li      t0, 0
    blt     t6, t0, .vbsub_next_col
    mul     t0, t6, s4
    add     t0, t0, s5
    slli    t0, t0, 3
    add     t0, s1, t0
    fld     ft0, 0(t0)
    addi    t5, t6, 1
.vbsub_k:
    bge     t5, s3, .vbsub_k_done
    mul     t1, t5, s3
    add     t1, t1, t6
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    mul     t2, t5, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fld     ft2, 0(t2)
    fnmsub.d ft0, ft1, ft2, ft0
    addi    t5, t5, 1
    j       .vbsub_k
.vbsub_k_done:
    mul     t1, t6, s3
    add     t1, t1, t6
    slli    t1, t1, 3
    add     t1, s0, t1
    fld     ft1, 0(t1)
    fdiv.d  ft0, ft0, ft1
    mul     t2, t6, s4
    add     t2, t2, s5
    slli    t2, t2, 3
    add     t2, s2, t2
    fsd     ft0, 0(t2)
    addi    t6, t6, -1
    j       .vbsub_row
.vbsub_next_col:
    addi    s5, s5, 1
    j       .vbsub_col
.vbsub_done:
    ld      ra,  0(sp)
    ld      s0,  8(sp)
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
    sd      ra,  0(sp)
    sd      s0,  8(sp)
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
    mv      a0, s0
    mv      a1, s4
    call    asm_cholesky
    mv      a0, s0
    mv      a1, s1
    mv      a2, s3
    mv      a3, s4
    mv      a4, s5
    call    asm_forward_sub
    mv      a0, s0
    mv      a1, s3
    mv      a2, s2
    mv      a3, s4
    mv      a4, s5
    call    asm_back_sub
    ld      ra,  0(sp)
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    addi    sp, sp, 48
    ret

# Read-only constants
    .section .rodata
    .align 3

.vec_half:
    .double 0.5

.vec_zero:
    .double 0.0

.vec_neg_one:
    .double -1.0

.vec_eps:
    .double 1e-15