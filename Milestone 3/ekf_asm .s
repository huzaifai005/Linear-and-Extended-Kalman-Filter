    .section .text
    .global asm_cordic_atan2
    .global asm_h_joint
    .global asm_h_full
    .global asm_jacobian_joint
    .global asm_jacobian_full
    .global asm_ekf_predict
    .global asm_ekf_update

asm_cordic_atan2:
    addi    sp, sp, -48
    sd      ra, 0(sp)
    fsd     fs0, 8(sp)
    fsd     fs1, 16(sp)
    fsd     fs2, 24(sp)
    fsd     fs3, 32(sp)
    fsd     fs4, 40(sp)

    # Load constants
    la      t3, .cordic_half
    fld     fa7, 0(t3)     

    la      t3, .cordic_zero
    fld     ft2, 0(t3)      

    la      t3, .cordic_pi
    fld     ft1, 0(t3)     

    la      t3, .cordic_eps
    fld     fs4, 0(t3)      

    # Special case: x ≈ 0 and y ≈ 0 (origin) -> return 0
    fmul.d  fs0, fa0, fa0   
    fmul.d  fs1, fa1, fa1   
    fadd.d  fs0, fs0, fs1   
    flt.d   t0, fs0, fs4    
    beqz    t0, .cordic_not_origin
    la      t3, .cordic_zero
    fld     fa0, 0(t3)
    j       .cordic_ret
.cordic_not_origin:

    # Special case: x ≈ 0 
    # If x is near zero: return PI/2 (y>0) or -PI/2 (y<0)
    flt.d   t0, fs1, fs4   
    beqz    t0, .cordic_x_not_zero
    la      t3, .cordic_pi_half
    fld     fa0, 0(t3)      
    flt.d   t0, fa0, ft2   
    # Re-do: check original y 
    fmv.d   fs1, fa0        
    la      t3, .cordic_pi_half
    fld     fa0, 0(t3)      
    flt.d   t0, fs1, ft2   
    beqz    t0, .cordic_ret 
    la      t3, .cordic_neg_pi_half
    fld     fa0, 0(t3)      
    j       .cordic_ret
.cordic_x_not_zero:

    # Quadrant handling 
    fmv.d   fs1, fa0        
    fmv.d   fs2, fa1        

    flt.d   t2, fs2, ft2    
    beqz    t2, .cordic_q1q4   

    # x < 0
    la      t3, .cordic_neg_one
    fld     ft2, 0(t3)
    fmul.d  fa0, fs1, ft2  
    fmul.d  fa1, fs2, ft2   
    j       .cordic_run

.cordic_q1q4:
    fmv.d   fa0, fs1        
    fmv.d   fa1, fs2        

.cordic_run:
    # CORDIC vectoring mode 
    # fa1 = working x
    # fa0 = working y
    # fa3 = accumulated angle = 0
    la      t3, .cordic_zero
    fld     fa3, 0(t3)    
    la      t3, .cordic_one
    fld     fa4, 0(t3)      

    la      t3, .cordic_table  
    li      t0, 0      
    li      t1, 52         

.cordic_iter:
    bge     t0, t1, .cordic_iter_done

    # Load CORDIC_TABLE[i] 
    slli    t4, t0, 3       
    add     t4, t3, t4
    fld     ft0, 0(t4)     

    # Decide rotation direction based on sign of y:
    #   y > 0: rotate clockwise (positive angle)
    #   y < 0: rotate counter-clockwise       
    la      t4, .cordic_zero
    fld     ft2, 0(t4)
    flt.d   t4, ft2, fa0  

    # Compute x * power_of_2 and y * power_of_2
    fmul.d  fa5, fa0, fa4   
    fmul.d  fa6, fa1, fa4   

    beqz    t4, .cordic_neg_rot 

    # Clockwise rotation (y > 0):
    #   x_new = x + y * 2^{-i}
    #   y_new = y - x * 2^{-i}
    #   angle += table[i]
    fadd.d  fs3, fa1, fa5   
    fsub.d  fs0, fa0, fa6   
    fadd.d  fa3, fa3, ft0   
    j       .cordic_update

.cordic_neg_rot:
    # Counter-clockwise rotation (y <= 0):
    #   x_new = x - y * 2^{-i}
    #   y_new = y + x * 2^{-i}
    #   angle -= table[i]
    fsub.d  fs3, fa1, fa5
    fadd.d  fs0, fa0, fa6
    fsub.d  fa3, fa3, ft0

.cordic_update:
    fmv.d   fa1, fs3       
    fmv.d   fa0, fs0        
    fmul.d  fa4, fa4, fa7  

    addi    t0, t0, 1
    j       .cordic_iter

.cordic_iter_done:
    # fa3 = accumulated angle

    # Quadrant correction 
    la      t3, .cordic_pi
    fld     ft1, 0(t3)
    la      t3, .cordic_zero
    fld     ft2, 0(t3)

    flt.d   t2, fs2, ft2    
    beqz    t2, .cordic_no_correction

    # x < 0
    flt.d   t4, fs1, ft2   
    bnez    t4, .cordic_Q3
  
    fadd.d  fa3, fa3, ft1
    j       .cordic_no_correction
.cordic_Q3:
   
    fsub.d  fa3, fa3, ft1

.cordic_no_correction:
    fmv.d   fa0, fa3        

.cordic_ret:
    fld     fs0, 8(sp)
    fld     fs1, 16(sp)
    fld     fs2, 24(sp)
    fld     fs3, 32(sp)
    fld     fs4, 40(sp)
    ld      ra, 0(sp)
    addi    sp, sp, 48
    ret


asm_h_joint:
    addi    sp, sp, -56
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    fsd     fs0, 24(sp)
    fsd     fs1, 32(sp)
    fsd     fs2, 40(sp)
    fsd     fs3, 48(sp)

    mv      s0, a0         
    mv      s1, a1          

    # Load px, py, pz
    fld     fs0,  0(s0)    
    fld     fs1, 32(s0)    
    fld     fs2, 64(s0)    

    # rho = sqrt(px^2 + py^2)
    fmul.d  fa0, fs0, fs0   
    fmul.d  fa1, fs1, fs1   
    fadd.d  fa0, fa0, fa1   
    call    asm_newton_sqrt  
    fmv.d   fs3, fa0       

    # r = sqrt(px^2 + py^2 + pz^2) = sqrt(rho^2 + pz^2)
    fmul.d  fa0, fs3, fs3   
    fmul.d  fa1, fs2, fs2  
    fadd.d  fa0, fa0, fa1  
    call    asm_newton_sqrt 

    # Store r in z[0]
    fsd     fa0, 0(s1)

    # Check r > epsilon to avoid division by zero
    la      t0, .ekf_eps
    fld     ft0, 0(t0)
    flt.d   t1, fa0, ft0    
    bnez    t1, .h_joint_degenerate

    # theta = atan2(py, px)
    fmv.d   fa0, fs1        
    fmv.d   fa1, fs0        
    call    asm_cordic_atan2
    fsd     fa0, 8(s1)      

    # phi = atan2(pz, rho)
    fmv.d   fa0, fs2        
    fmv.d   fa1, fs3        
    call    asm_cordic_atan2
    fsd     fa0, 16(s1)     
    j       .h_joint_done

.h_joint_degenerate:
    la      t0, .ekf_zero
    fld     ft0, 0(t0)
    fsd     ft0, 8(s1)
    fsd     ft0, 16(s1)

.h_joint_done:
    fld     fs0, 24(sp)
    fld     fs1, 32(sp)
    fld     fs2, 40(sp)
    fld     fs3, 48(sp)
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      ra, 0(sp)
    addi    sp, sp, 56
    ret

asm_h_full:
    addi    sp, sp, -32
    sd      ra, 0(sp)
    sd      s0, 8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)

    mv      s0, a0          
    mv      s1, a1          
    li      s2, 0           

.h_full_loop:
    li      t0, 23
    bge     s2, t0, .h_full_done

    # Joint j state starts at offset j*12*8 = j*96 bytes
    li      t1, 96
    mul     t1, s2, t1
    add     a0, s0, t1    

    # Output for joint j starts at offset j*3*8 = j*24 bytes
    li      t2, 24
    mul     t2, s2, t2
    add     a1, s1, t2      

    call    asm_h_joint

    addi    s2, s2, 1
    j       .h_full_loop

.h_full_done:
    ld      s0, 8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      ra, 0(sp)
    addi    sp, sp, 32
    ret

asm_jacobian_joint:
    addi    sp, sp, -80
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    fsd     fs0, 24(sp)
    fsd     fs1, 32(sp)
    fsd     fs2, 40(sp)
    fsd     fs3, 48(sp)
    fsd     fs4, 56(sp)
    fsd     fs5, 64(sp)

    mv      s0, a0          # joint state
    mv      s1, a1          # output Jacobian

    # Zero the output matrix
    la      t0, .ekf_zero
    fld     ft0, 0(t0)
    li      t1, 0
    li      t2, 36          
.jac_zero:
    bge     t1, t2, .jac_zero_done
    slli    t3, t1, 3
    add     t3, s1, t3
    fsd     ft0, 0(t3)
    addi    t1, t1, 1
    j       .jac_zero
.jac_zero_done:

    # Load px, py, pz
    fld     fs0,  0(s0)    
    fld     fs1, 32(s0)    
    fld     fs2, 64(s0)    

    # rho2 = px^2 + py^2
    fmul.d  fs3, fs0, fs0
    fmul.d  fs4, fs1, fs1
    fadd.d  fs3, fs3, fs4   

    # r2 = rho2 + pz^2
    fmul.d  fs4, fs2, fs2
    fadd.d  fs4, fs4, fs3   

    # rho = sqrt(rho2)
    fmv.d   fa0, fs3
    call    asm_newton_sqrt
    fmv.d   fs5, fa0        

    # r = sqrt(r2)
    fmv.d   fa0, fs4
    call    asm_newton_sqrt


    # Check singularities
    la      t0, .ekf_eps
    fld     ft1, 0(t0)
    flt.d   t1, fa0, ft1    
    bnez    t1, .jac_done  
    flt.d   t2, fs5, ft1    

    # Row 0: dr/dx = px/r, py/r, pz/r
    # H[0][0] = px / r
    fdiv.d  ft2, fs0, fa0
    fsd     ft2, 0(s1)     

    # H[0][4] = py / r 
    fdiv.d  ft2, fs1, fa0
    fsd     ft2, 32(s1)    

    # H[0][8] = pz / r
    fdiv.d  ft2, fs2, fa0
    fsd     ft2, 64(s1)    

    bnez    t2, .jac_done

    # Row 1: dtheta/dx
    # H[1][0] = -py / rho2
    fdiv.d  ft2, fs1, fs3   
    la      t0, .ekf_neg_one
    fld     ft3, 0(t0)
    fmul.d  ft2, ft2, ft3   
    fsd     ft2, 96(s1)  

    # H[1][4] = px / rho2
    fdiv.d  ft2, fs0, fs3
    fsd     ft2, 128(s1)    

    # Row 2: dphi/dx
    # denominator_xz = rho * r2
    fmul.d  ft4, fs5, fs4   

    # H[2][0] = -(px * pz) / (rho * r2)
    fmul.d  ft2, fs0, fs2  
    fdiv.d  ft2, ft2, ft4  
    fmul.d  ft2, ft2, ft3  
    fsd     ft2, 192(s1)   

    # H[2][4] = -(py * pz) / (rho * r2)
    fmul.d  ft2, fs1, fs2
    fdiv.d  ft2, ft2, ft4
    fmul.d  ft2, ft2, ft3
    fsd     ft2, 224(s1)   

    # H[2][8] = rho / r2
    fdiv.d  ft2, fs5, fs4
    fsd     ft2, 256(s1)    

.jac_done:
    fld     fs0, 24(sp)
    fld     fs1, 32(sp)
    fld     fs2, 40(sp)
    fld     fs3, 48(sp)
    fld     fs4, 56(sp)
    fld     fs5, 64(sp)
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      ra,  0(sp)
    addi    sp, sp, 80
    ret

asm_jacobian_full:
    addi    sp, sp, -48
    sd      ra,  0(sp)
    sd      s0,  8(sp)
    sd      s1, 16(sp)
    sd      s2, 24(sp)
    sd      s3, 32(sp)

    mv      s0, a0          
    mv      s1, a1          
    li      s2, 0           

.jac_full_loop:
    li      t0, 23
    bge     s2, t0, .jac_full_done

    # Pointer to joint j's state: s0 + j*12*8 = s0 + j*96
    li      t1, 96
    mul     t1, s2, t1
    add     a0, s0, t1

    li      t2, 840
    mul     t2, s2, t2
    slli    t2, t2, 3
    add     a1, s1, t2

    mv      s3, a2        
    mv      a1, s3          
    mv      a2, a1         
    call    asm_jacobian_joint

    # Copy temp_3x12 into Hk at the right block with row stride 276
    li      t3, 0          
.copy_row:
    li      t4, 3
    bge     t3, t4, .copy_done

    li      t5, 0          
.copy_col:
    li      t6, 12
    bge     t5, t6, .next_row

    mul     t0, t3, t6      
    add     t0, t0, t5      
    slli    t0, t0, 3
    add     t0, s3, t0
    fld     ft0, 0(t0)

    # Dest: Hk[(j*3+r)*276 + (j*12+c)]
    li      t1, 3
    mul     t1, s2, t1      
    add     t1, t1, t3      
    li      t2, 276
    mul     t1, t1, t2      
    li      t2, 12
    mul     t2, s2, t2      
    add     t1, t1, t2      
    add     t1, t1, t5      
    slli    t1, t1, 3
    add     t1, s1, t1
    fsd     ft0, 0(t1)

    addi    t5, t5, 1
    j       .copy_col
.next_row:
    addi    t3, t3, 1
    j       .copy_row
.copy_done:
    addi    s2, s2, 1
    j       .jac_full_loop

.jac_full_done:
    ld      s0,  8(sp)
    ld      s1, 16(sp)
    ld      s2, 24(sp)
    ld      s3, 32(sp)
    ld      ra,  0(sp)
    addi    sp, sp, 48
    ret

    .extern asm_lkf_predict
asm_ekf_predict:
    la      t1, asm_lkf_predict
    jr      t1                 

asm_ekf_update:
    addi    sp, sp, -144
    sd      ra,    0(sp)
    sd      s0,    8(sp)
    sd      s1,   16(sp)
    sd      s2,   24(sp)
    sd      s3,   32(sp)
    sd      s4,   40(sp)
    sd      s5,   48(sp)
    sd      s6,   56(sp)
    sd      s7,   64(sp)
    sd      s8,   72(sp)
    sd      s9,   80(sp)
    sd      s10,  88(sp)
    sd      s11,  96(sp)
    fsd     fs8,  104(sp)
    fsd     fs9,  112(sp)
    fsd     fs10, 120(sp)

    mv      s0, a0          
    mv      s1, a1          
    mv      s2, a2         
    mv      s4, a4          
    mv      s5, a5          
    mv      s6, a6          
    fmv.d   fs8, fa0        

    ld      s8, 144(sp)     
    ld      s9, 152(sp)     

    # Workspace pointers
    mv      s10, s6

    mul     t0, s8, s8
    slli    t0, t0, 3
    add     s11, s6, t0    

    # Build adaptive R_ekf 
    mul     t0, s8, s8      
    li      t1, 0
    la      t2, .ekf_zero
    fld     ft0, 0(t2)
.zero_R:
    bge     t1, t0, .zero_R_done
    slli    t2, t1, 3
    add     t2, s10, t2
    fsd     ft0, 0(t2)
    addi    t1, t1, 1
    j       .zero_R
.zero_R_done:

    fmul.d  fs9, fs8, fs8   

    la      t0, .ekf_min_rho2
    fld     fs10, 0(t0)     

    li      t3, 0           
.build_R:
    li      t4, 23
    bge     t3, t4, .R_done

    # Load px, py, pz for joint j
    li      t0, 96
    mul     t0, t3, t0
    add     t0, s0, t0     
    fld     fa0,  0(t0)     
    fld     fa1, 32(t0)    
    fld     fa2, 64(t0)     

    # rho2 = max(px^2 + py^2, min_rho2)
    fmul.d  ft0, fa0, fa0
    fmul.d  ft1, fa1, fa1
    fadd.d  ft0, ft0, ft1   
    flt.d   t5, ft0, fs10   
    beqz    t5, .rho2_ok
    fmv.d   ft0, fs10       
.rho2_ok:
    # r2 = rho2 + pz^2
    fmul.d  ft2, fa2, fa2
    fadd.d  ft1, ft0, ft2   
    flt.d   t5, ft1, fs10
    beqz    t5, .r2_ok
    fmv.d   ft1, fs10
.r2_ok:
    # Base offset in R
    li      t5, 3
    mul     t5, t3, t5  
    # R_range: R[j*3][j*3] = sigma^2
    mul     t6, t5, s8
    add     t6, t6, t5    
    slli    t6, t6, 3
    add     t6, s10, t6
    fsd     fs9, 0(t6)    

    # R_azimuth: R[j*3+1][j*3+1] = sigma^2 / rho2
    fdiv.d  ft3, fs9, ft0
    addi    t5, t5, 1    
    mul     t6, t5, s8
    add     t6, t6, t5
    slli    t6, t6, 3
    add     t6, s10, t6
    fsd     ft3, 0(t6)

    # R_elevation: R[j*3+2][j*3+2] = sigma^2 / r2
    fdiv.d  ft3, fs9, ft1
    addi    t5, t5, 1    
    mul     t6, t5, s8
    add     t6, t6, t5
    slli    t6, t6, 3
    add     t6, s10, t6
    fsd     ft3, 0(t6)

    addi    t3, t3, 1
    j       .build_R
.R_done:

    # Step 2: Hk = jacobian_full(x_pred)
    mul     t0, s8, s7     
    li      t1, 0
    la      t2, .ekf_zero
    fld     ft0, 0(t2)
.zero_Hk:
    bge     t1, t0, .zero_Hk_done
    slli    t2, t1, 3
    add     t2, s11, t2
    fsd     ft0, 0(t2)
    addi    t1, t1, 1
    j       .zero_Hk
.zero_Hk_done:

    mv      a0, s0          
    mv      a1, s11         
    # Push temp_3x12 on stack for jacobian_full
    addi    sp, sp, -8
    sd      s9, 0(sp)
    call    asm_jacobian_full
    addi    sp, sp, 8

    # Step 3: z_hat = h_full(x_pred)
    # Store z_hat after Hk in workbuf
    mul     t0, s8, s7      
    slli    t0, t0, 3
    add     s3, s11, t0     

    mv      a0, s0          
    mv      a1, s3         
    call    asm_h_full

    # Step 4: y = z_sph - z_hat; wrap angles
    # y stored after z_hat
    slli    t0, s8, 3       
    add     t1, s3, t0      

    mv      a0, s2          
    mv      a1, s3          
    mv      a2, t1          
    mv      a3, s8          
    li      a4, 1
    call    asm_mat_sub

    # Angle wrapping: for each joint j, wrap y[j*3+1] and y[j*3+2]
    li      t3, 0
.wrap_loop:
    li      t4, 23
    bge     t3, t4, .wrap_done

    # azimuth index = j*3 + 1, elevation = j*3 + 2
    li      t5, 3
    mul     t5, t3, t5    
    addi    t6, t5, 1      

    # Load y[j*3+1]
    slli    t0, t6, 3
    add     t0, t1, t0
    fld     fa0, 0(t0)
    call    asm_wrap_angle
    fsd     fa0, 0(t0)

    # Load y[j*3+2] 
    addi    t6, t6, 1
    slli    t0, t6, 3
    add     t0, t1, t0
    fld     fa0, 0(t0)
    call    asm_wrap_angle
    fsd     fa0, 0(t0)

    addi    t3, t3, 1
    j       .wrap_loop
.wrap_done:

    fld     fs8,  104(sp)
    fld     fs9,  112(sp)
    fld     fs10, 120(sp)
    ld      s0,    8(sp)
    ld      s1,   16(sp)
    ld      s2,   24(sp)
    ld      s3,   32(sp)
    ld      s4,   40(sp)
    ld      s5,   48(sp)
    ld      s6,   56(sp)
    ld      s7,   64(sp)
    ld      s8,   72(sp)
    ld      s9,   80(sp)
    ld      s10,  88(sp)
    ld      s11,  96(sp)
    ld      ra,    0(sp)
    addi    sp, sp, 144
    ret

    .global asm_wrap_angle
asm_wrap_angle:
    la      t0, .cordic_pi
    fld     ft0, 0(t0)          
    la      t0, .cordic_two_pi
    fld     ft1, 0(t0)          
.wrap_pos:
    flt.d   t1, ft0, fa0        
    beqz    t1, .wrap_neg
    fsub.d  fa0, fa0, ft1       
    j       .wrap_pos
.wrap_neg:
    la      t0, .cordic_neg_pi
    fld     ft2, 0(t0)
    flt.d   t1, fa0, ft2       
    beqz    t1, .wrap_done_angle
    fadd.d  fa0, fa0, ft1     
    j       .wrap_neg
.wrap_done_angle:
    ret

# Read-only constants for ekf_asm.s

    .section .rodata
    .align 3

.ekf_zero:      .double 0.0
.ekf_eps:       .double 1e-9
.ekf_neg_one:   .double -1.0
.ekf_min_rho2:  .double 0.01

# CORDIC constants
.cordic_zero:       .double 0.0
.cordic_one:        .double 1.0
.cordic_half:       .double 0.5
.cordic_neg_one:    .double -1.0
.cordic_eps:        .double 1e-9
.cordic_pi:         .double 3.141592653589793
.cordic_neg_pi:     .double -3.141592653589793
.cordic_pi_half:    .double 1.5707963267948966
.cordic_neg_pi_half:.double -1.5707963267948966
.cordic_two_pi:     .double 6.283185307179586

# CORDIC lookup table: atan(2^{-i}) for i = 0..51
# Each entry is a double (8 bytes).
.cordic_table:
    .double  0.7853981633974483        # i=0:  atan(2^0)  = pi/4
    .double  0.4636476090008257        # i=1:  atan(0.5)
    .double  0.24497866312686414       # i=2:  atan(0.25)
    .double  0.12435499454676144       # i=3
    .double  0.06241880999595735       # i=4
    .double  0.031239833430268277      # i=5
    .double  0.015623728620476831      # i=6
    .double  0.007812341060101111      # i=7
    .double  0.003906230131966972      # i=8
    .double  0.0019531225164788188     # i=9
    .double  0.0009765621895593195     # i=10
    .double  0.0004882812111948983     # i=11
    .double  0.00024414062014936177    # i=12
    .double  0.00012207031189367021    # i=13
    .double  0.000061035156174208773   # i=14
    .double  0.000030517578115526096   # i=15
    .double  0.000015258789062277982   # i=16
    .double  0.0000076293945311369800  # i=17
    .double  0.0000038146972656429850  # i=18
    .double  0.0000019073486329385070  # i=19
    .double  0.00000095367431640625    # i=20
    .double  0.000000476837158203125   # i=21
    .double  0.0000002384185791015625  # i=22
    .double  0.00000011920928955078125 # i=23
    # i=24..51: atan(2^{-i}) == 2^{-i} to machine precision
    .double  5.9604644775390625e-08    # i=24: 2^{-24}
    .double  2.9802322387695312e-08    # i=25: 2^{-25}
    .double  1.4901161193847656e-08    # i=26: 2^{-26}
    .double  7.450580596923828e-09     # i=27: 2^{-27}
    .double  3.725290298461914e-09     # i=28: 2^{-28}
    .double  1.862645149230957e-09     # i=29: 2^{-29}
    .double  9.313225746154785e-10     # i=30: 2^{-30}
    .double  4.656612873077393e-10     # i=31: 2^{-31}
    .double  2.3283064365386963e-10    # i=32: 2^{-32}
    .double  1.1641532182693481e-10    # i=33: 2^{-33}
    .double  5.820766091346741e-11     # i=34: 2^{-34}
    .double  2.9103830456733704e-11    # i=35: 2^{-35}
    .double  1.4551915228366852e-11    # i=36: 2^{-36}
    .double  7.275957614183426e-12     # i=37: 2^{-37}
    .double  3.637978807091713e-12     # i=38: 2^{-38}
    .double  1.8189894035458565e-12    # i=39: 2^{-39}
    .double  9.094947017729282e-13     # i=40: 2^{-40}
    .double  4.547473508864641e-13     # i=41: 2^{-41}
    .double  2.2737367544323206e-13    # i=42: 2^{-42}
    .double  1.1368683772161603e-13    # i=43: 2^{-43}
    .double  5.684341886080802e-14     # i=44: 2^{-44}
    .double  2.842170943040401e-14     # i=45: 2^{-45}
    .double  1.4210854715202004e-14    # i=46: 2^{-46}
    .double  7.105427357601002e-15     # i=47: 2^{-47}
    .double  3.552713678800501e-15     # i=48: 2^{-48}
    .double  1.7763568394002505e-15    # i=49: 2^{-49}
    .double  8.881784197001252e-16     # i=50: 2^{-50}
    .double  4.440892098500626e-16     # i=51: 2^{-51} ~= machine epsilon