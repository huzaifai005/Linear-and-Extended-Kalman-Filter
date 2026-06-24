    .section .text


    .global asm_vec_zero_matrix
asm_vec_zero_matrix:
    mv      t0, a1              # total elements remaining
    mv      t1, a0              # pointer
    fmv.d.x ft0, zero           # ft0 = 0.0
.vzero_loop:
    beqz    t0, .vzero_done
    vsetvli t2, t0, e64, m4, ta, ma
    vfmv.v.f v0, ft0            # v0..v3 = all 0.0
    vse64.v  v0, (t1)
    slli    t3, t2, 3
    add     t1, t1, t3
    sub     t0, t0, t2
    j       .vzero_loop
.vzero_done:
    ret


    .global asm_vec_mat_copy
asm_vec_mat_copy:
    mv      t0, a2
    mv      t1, a0              # dst
    mv      t2, a1              # src
.vcopy_loop:
    beqz    t0, .vcopy_done
    vsetvli t3, t0, e64, m4, ta, ma
    vle64.v v0, (t2)
    vse64.v v0, (t1)
    slli    t4, t3, 3
    add     t1, t1, t4
    add     t2, t2, t4
    sub     t0, t0, t3
    j       .vcopy_loop
.vcopy_done:
    ret


    .global asm_vec_negate
asm_vec_negate:
    mv      t0, a1
    mv      t1, a0
    la      t2, .ekfv_neg_one
    fld     ft0, 0(t2)          # ft0 = -1.0
.vneg_loop:
    beqz    t0, .vneg_done
    vsetvli t2, t0, e64, m4, ta, ma
    vle64.v v0, (t1)
    vfmul.vf v0, v0, ft0        # v0 = -v0
    vse64.v v0, (t1)
    slli    t3, t2, 3
    add     t1, t1, t3
    sub     t0, t0, t2
    j       .vneg_loop
.vneg_done:
    ret


    .global asm_vec_add_identity
asm_vec_add_identity:
    # stride = (N+1)*8
    addi    t0, a1, 1           # N+1
    slli    t0, t0, 3           # stride = (N+1)*8
    mv      t1, a0              # pointer to A[0][0]
    mv      t2, a1              # remaining = N
    la      t3, .ekfv_one
    fld     ft0, 0(t3)          # ft0 = 1.0
.vdiag_loop:
    beqz    t2, .vdiag_done
    fld     ft1, 0(t1)
    fadd.d  ft1, ft1, ft0
    fsd     ft1, 0(t1)
    add     t1, t1, t0          # advance by stride
    addi    t2, t2, -1
    j       .vdiag_loop
.vdiag_done:
    ret

# Read-only constants
    .section .rodata
    .align 3

.ekfv_neg_one:
    .double -1.0

.ekfv_one:
    .double 1.0

.ekfv_zero:
    .double 0.0