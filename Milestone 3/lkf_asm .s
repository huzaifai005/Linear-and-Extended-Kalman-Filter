    .section .text
    .global asm_lkf_predict
    .global asm_lkf_update

asm_lkf_predict:
    # Save callee-saved registers
    addi    sp, sp, -80
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)
    sd      s4, 40(sp)
    sd      s5, 48(sp)
    sd      s6, 56(sp)
    sd      s7, 64(sp)

    # Move arguments to saved registers
    mv      s0, a0          
    mv      s1, a1          
    mv      s2, a2          
    mv      s3, a3          
    mv      s4, a4          
    mv      s5, a5          
    mv      s6, a6          
    mv      s7, a7          

    # Step 1: x_pred = F * x 
    # mat_mul(F, x, x_pred, N, N, 1)
    mv      a0, s2         
    mv      a1, s0         
    mv      a2, s4        
    mv      a3, s7          
    mv      a4, s7         
    li      a5, 1           
    call    asm_mat_mul

    # Step 2: work1 = F * P 
    # mat_mul(F, P, work1, N, N, N)
    mv      a0, s2          
    mv      a1, s1          
    mv      a2, s6          
    mv      a3, s7          
    mv      a4, s7          
    mv      a5, s7          
    call    asm_mat_mul

    # Step 3: P_pred = work1 * F^T 
    mv      a0, s2          
    mv      a1, s5         
    mv      a2, s7          
    mv      a3, s7          
    call    asm_mat_transpose   

    ld      t0, 80(sp)    

    # Compute: work2 = work1 * F^T
    # F^T is currently stored in P_pred (s5)
    mv      a0, s6          
    mv      a1, s5         
    mv      a2, t0          
    mv      a3, s7          
    mv      a4, s7          
    mv      a5, s7          
    call    asm_mat_mul
    # work2 now holds F * P * F^T

    # Step 4: P_pred = (F*P*F^T) + Q
    mv      a0, t0          
    mv      a1, s3         
    mv      a2, s5          
    mv      a3, s7        
    mv      a4, s7          
    call    asm_mat_add

    # Step 5: symmetrize P_pred
    mv      a0, s5         
    mv      a1, s7         
    call    asm_symmetrize

    # Restore and return
    ld      ra,  0(sp)
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      s4, 40(sp)
    ld      s5, 48(sp)
    ld      s6, 56(sp)
    ld      s7, 64(sp)
    addi    sp, sp, 80
    ret

asm_lkf_update:
    addi    sp, sp, -112
    sd      ra,   0(sp)
    sd      s0,   8(sp)
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
    sd      s11, 96(sp)

    mv      s0, a0          
    mv      s1, a1          
    mv      s2, a2         
    mv      s3, a3         
    mv      s4, a4          
    mv      s5, a5         
    mv      s6, a6         
    mv      s7, a7         

    # Load stack-passed N and M
    ld      s8, 112(sp)     
    ld      s9, 120(sp)     

    # Compute workspace sub-pointers 
    # S_copy  at workbuf + 0
    mv      s10, s7         

    # PHt     at workbuf + M*M*8
    mul     t0, s9, s9
    slli    t0, t0, 3
    add     s11, s7, t0     

    # Kt      at workbuf + (M*M + N*M)*8
    mul     t1, s8, s9
    slli    t1, t1, 3
    add     t2, s11, t1     

    # Hxp     at workbuf + (M*M + 2*N*M)*8
    add     t3, t2, t1      

    # IKH     at workbuf + (M*M + 2*N*M + M)*8
    slli    t4, s9, 3       
    add     t4, t3, t4      

    # tmp1    at IKH + N*N*8
    mul     t5, s8, s8
    slli    t5, t5, 3
    add     t5, t4, t5      

    # tmp2    at tmp1 + N*N*8
    mul     t6, s8, s8
    slli    t6, t6, 3
    add     t6, t5, t6      

    
    addi    sp, sp, -56     
    sd      s11,  0(sp)     
    sd      t2,   8(sp)     
    sd      t3,  16(sp)     
    sd      t4,  24(sp)     
    sd      t5,  32(sp)     
    sd      t6,  40(sp)    
    sd      s10, 48(sp)     

    
    ld      t4, 24(sp)      
    mv      a0, s2         
    mv      a1, t4          
    mv      a2, s9          
    mv      a3, s8          
    call    asm_mat_transpose

    ld      s11, 0(sp)      
    mv      a0, s1          
    mv      a1, t4          
    mv      a2, s11         
    mv      a3, s8          
    mv      a4, s8          
    mv      a5, s9          
    call    asm_mat_mul    

    # Step 2: S = H * PHt + R
    ld      s10, 48(sp)     
    mv      a0, s2         
    mv      a1, s11         
    mv      a2, s10         
    mv      a3, s9          
    mv      a4, s8          
    mv      a5, s9          
    call    asm_mat_mul     

    mv      a0, s10         
    mv      a1, s3          
    mv      a2, s10         
    mv      a3, s9          
    mv      a4, s9          
    call    asm_mat_add     

    mv      a0, s10
    mv      a1, s9
    call    asm_symmetrize  

    # Step 3: Kt = S^{-1} * PHt^T  (Cholesky solve)
    ld      t2, 8(sp)      
    mv      a0, s11         
    mv      a1, t2          
    mv      a2, s8          
    mv      a3, s9         
    call    asm_mat_transpose  

    # Allocate Ychol temp: use tmp2 as workspace
    ld      t6, 40(sp)      
    mv      a0, s10        
    mv      a1, t2          
    mv      a2, t2         
                            
    mv      a3, t6        
    mv      a4, s9          
    mv      a5, s8          
    call    asm_cholesky_solve
    # t2 (Kt) now holds K^T (MxN), so K = Kt^T is NxM

    # Step 4: K = Kt^T  (NxM) — store into tmp1
    ld      t5, 32(sp)      
    mv      a0, t2         
    mv      a1, t5         
    mv      a2, s9          
    mv      a3, s8         
    call    asm_mat_transpose  

    # Step 5: y = z - H * x_pred
    ld      t3, 16(sp)     
    mv      a0, s2          
    mv      a1, s0          
    mv      a2, t3          
    mv      a3, s9          
    mv      a4, s8          
    li      a5, 1
    call    asm_mat_mul     

    # y = z - Hxp  (stored back in Hxp to save memory)
    mv      a0, s4          
    mv      a1, t3          
    mv      a2, t3          
    mv      a3, s9          
    li      a4, 1
    call    asm_mat_sub     

    # Step 6: x_upd = x_pred + K * y
    # K is at tmp1 (NxM), y is at Hxp (Mx1)
    mv      a0, t5         
    mv      a1, t3          
    mv      a2, s5          
    mv      a3, s8          
    mv      a4, s9          
    li      a5, 1
    call    asm_mat_mul     

    mv      a0, s0          
    mv      a1, s5          
    mv      a2, s5          
    mv      a3, s8
    li      a4, 1
    call    asm_mat_add     

    # Step 7: IKH = I - K * H
    ld      t4, 24(sp)      
    mv      a0, t5          
    mv      a1, s2          
    mv      a2, t4          
    mv      a3, s8          
    mv      a4, s9          
    mv      a5, s8          
    call    asm_mat_mul     

    # IKH = I - K*H  (negate all elements then add I on diagonal)
    mul     t0, s8, s8
    li      t1, 0
    # negate: IKH[k] = -IKH[k]
.negate_KH:
    bge     t1, t0, .negate_done
    slli    t2, t1, 3
    add     t2, t4, t2
    fld     ft0, 0(t2)
    la      t3, .lkf_neg_one
    fld     ft1, 0(t3)
    fmul.d  ft0, ft0, ft1
    fsd     ft0, 0(t2)
    addi    t1, t1, 1
    j       .negate_KH
.negate_done:
    # Add 1.0 to diagonal
    la      t3, .lkf_one
    fld     ft1, 0(t3)
    li      t1, 0
.add_diag:
    bge     t1, s8, .add_diag_done
    mul     t2, t1, s8
    add     t2, t2, t1
    slli    t2, t2, 3
    add     t2, t4, t2
    fld     ft0, 0(t2)
    fadd.d  ft0, ft0, ft1
    fsd     ft0, 0(t2)
    addi    t1, t1, 1
    j       .add_diag
.add_diag_done:

    # Step 8: P_upd = IKH * P_pred * IKH^T + K*R*K^T (Joseph form)
    # tmp2 = IKH * P_pred
    ld      t6, 40(sp)      
    mv      a0, t4          
    mv      a1, s1         
    mv      a2, t6         
    mv      a3, s8
    mv      a4, s8
    mv      a5, s8
    call    asm_mat_mul

    # IKH^T stored in tmp1 (reusing, K no longer needed)
    mv      a0, t4          
    mv      a1, t5          
    mv      a2, s8
    mv      a3, s8
    call    asm_mat_transpose

    # P_upd = tmp2 * IKH^T
    mv      a0, t6         
    mv      a1, t5          
    mv      a2, s6         
    mv      a3, s8
    mv      a4, s8
    mv      a5, s8
    call    asm_mat_mul

    # Now add K*R*K^T:
    # tmp2 = K * R  (NxM * MxM = NxM)
    ld      t2, 8(sp)       
    ld      t5, 32(sp)      
    mv      a0, t5          
    mv      a1, s3          
    mv      a2, t2          
    mv      a3, s8
    mv      a4, s9
    mv      a5, s9
    call    asm_mat_mul

    # K^T (MxN) — store in tmp1
    mv      a0, t5          
    mv      a1, t5          
    # Use IKH space temporarily: we already applied IKH, can overwrite it.
    mv      a0, t5        
    mv      a1, t4          
    mv      a2, s8         
    mv      a3, s9          
    call    asm_mat_transpose

    # tmp2_2 = (K*R) * K^T (NxM * MxN = NxN) — store in tmp1 space
    mv      a0, t2          
    mv      a1, t4          
    mv      a2, t5          
    mv      a3, s8
    mv      a4, s9
    mv      a5, s8
    call    asm_mat_mul

    # P_upd += K*R*K^T
    mv      a0, s6          
    mv      a1, t5          
    mv      a2, s6         
    mv      a3, s8
    mv      a4, s8
    call    asm_mat_add

    # Symmetrize P_upd
    mv      a0, s6
    mv      a1, s8
    call    asm_symmetrize

    # Epilogue
    addi    sp, sp, 56      
    ld      ra,   0(sp)
    ld      s0,   8(sp)
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
    ld      s11, 96(sp)
    addi    sp, sp, 112
    ret

# Read-only constants for lkf_asm.s
    .section .rodata
    .align 3
.lkf_neg_one:
    .double -1.0
.lkf_one:
    .double  1.0