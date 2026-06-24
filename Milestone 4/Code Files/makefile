
# Toolchain
CC      = riscv64-linux-gnu-gcc
AS      = riscv64-linux-gnu-gcc
LD      = riscv64-linux-gnu-gcc

# Runners
RUNNER      = qemu-riscv64 -cpu rv64,v=true,vlen=256,vext_spec=v1.0
RUNNER_VEC  = qemu-riscv64 -cpu rv64,v=true,vlen=128,vext_spec=v1.0

# Scalar flags (M3)
CFLAGS      = -march=rv64imafd -mabi=lp64d -O0 -g -Wall
ASFLAGS     = -march=rv64imafd -mabi=lp64d -g
LDFLAGS     = -march=rv64imafd -mabi=lp64d -static -lm

# Vector flags (M4) 
CFLAGS_VEC  = -march=rv64gcv -mabi=lp64d -O0 -g -Wall
ASFLAGS_VEC = -march=rv64gcv -mabi=lp64d -g
LDFLAGS_VEC = -march=rv64gcv -mabi=lp64d -static -lm

# Source files
ASM_COMMON_SCALAR = kalman_common.s
ASM_HELPERS_SCALAR = kalman_helpers.s
ASM_LKF      = lkf_asm.s
ASM_EKF      = ekf_asm.s
ASM_VEC_LKF  = lkf_vector.s
ASM_VEC_EKF  = ekf_vector.s
DRV_LKF      = driver_lkf.c
DRV_EKF      = driver_ekf.c

# Object files for scalar M3
OBJ_COMMON_SCALAR = kalman_common.o
OBJ_HELPERS_SCALAR = kalman_helpers.o
OBJ_LKF_ASM   = lkf_asm.o
OBJ_EKF_ASM   = ekf_asm.o
OBJ_DRV_LKF_SCALAR = driver_lkf_scalar.o
OBJ_DRV_EKF_SCALAR = driver_ekf_scalar.o

# Object files for vector M4 
OBJ_VEC_LKF   = lkf_vector.o
OBJ_VEC_EKF   = ekf_vector.o
OBJ_LKF_ASM_VEC = lkf_asm.o       
OBJ_EKF_ASM_VEC = ekf_asm.o        
OBJ_DRV_LKF_VEC = driver_lkf_vec.o
OBJ_DRV_EKF_VEC = driver_ekf_vec.o

# Data files 
NOISY_CSV   = noisy_data.csv
M2_LKF_OUT  = lkf_output.csv
M2_EKF_OUT  = ekf_output_final.csv

# Output files 
LKF_OUT     = lkf_asm_output.csv
EKF_OUT     = ekf_asm_output.csv
LKF_VEC_OUT = lkf_vec_output.csv
EKF_VEC_OUT = ekf_vec_output.csv

# Top-level targets

.PHONY: all all_vec all_m4 clean \
        run_lkf run_ekf run_lkf_vec run_ekf_vec \
        verify verify_vec help

all: lkf_asm ekf_asm
	@echo ""
	@echo "  Scalar build complete (Milestone 3)."
	@echo "  make run_lkf      — run scalar LKF"
	@echo "  make run_ekf      — run scalar EKF"

all_vec: lkf_vec ekf_vec
	@echo ""
	@echo "  Vector build complete (Milestone 4)."
	@echo "  make run_lkf_vec  — run vector LKF"
	@echo "  make run_ekf_vec  — run vector EKF"

all_m4: lkf_asm ekf_asm lkf_vec ekf_vec
	@echo ""
	@echo "  Full build complete (M3 scalar + M4 vector)."

# MILESTONE 3 — Scalar build rules

kalman_common.o: $(ASM_COMMON_SCALAR)
	$(AS) $(ASFLAGS) -c $< -o $@

kalman_helpers.o: $(ASM_HELPERS_SCALAR)
	$(AS) $(ASFLAGS) -c $< -o $@

lkf_asm.o: $(ASM_LKF)
	$(AS) $(ASFLAGS) -c $< -o $@

ekf_asm.o: $(ASM_EKF)
	$(AS) $(ASFLAGS) -c $< -o $@

driver_lkf_scalar.o: $(DRV_LKF)
	$(CC) $(CFLAGS) -c $< -o $@

driver_ekf_scalar.o: $(DRV_EKF)
	$(CC) $(CFLAGS) -c $< -o $@

lkf_asm: kalman_common.o kalman_helpers.o lkf_asm.o driver_lkf_scalar.o
	$(LD) $(LDFLAGS) $^ -o $@ -lm
	@echo "  Linked: lkf_asm (scalar)"

ekf_asm: kalman_common.o kalman_helpers.o ekf_asm.o lkf_asm.o driver_ekf_scalar.o
	$(LD) $(LDFLAGS) $^ -o $@ -lm
	@echo "  Linked: ekf_asm (scalar)"

# MILESTONE 4 — Vector build rules

lkf_vector.o: $(ASM_VEC_LKF)
	$(AS) $(ASFLAGS_VEC) -c $< -o $@

ekf_vector.o: $(ASM_VEC_EKF)
	$(AS) $(ASFLAGS_VEC) -c $< -o $@

driver_lkf_vec.o: $(DRV_LKF)
	$(CC) $(CFLAGS_VEC) -c $< -o $@

driver_ekf_vec.o: $(DRV_EKF)
	$(CC) $(CFLAGS_VEC) -c $< -o $@

lkf_vec: lkf_vector.o lkf_asm.o driver_lkf_vec.o
	$(LD) $(LDFLAGS_VEC) $^ -o $@ -lm
	@echo "  Linked: lkf_vec (vector)"


ekf_vec: ekf_vector.o lkf_vector.o ekf_asm.o lkf_asm.o driver_ekf_vec.o
	$(LD) $(LDFLAGS_VEC) $^ -o $@ -lm
	@echo "  Linked: ekf_vec (vector)"

# Run targets — Milestone 3
run_lkf: lkf_asm
	@echo "--- Running Scalar LKF (M3) ---"
	$(RUNNER) ./lkf_asm $(NOISY_CSV) $(LKF_OUT)

run_ekf: ekf_asm
	@echo "--- Running Scalar EKF (M3) ---"
	$(RUNNER) ./ekf_asm $(NOISY_CSV) $(EKF_OUT)

# Run targets — Milestone 4


run_lkf_vec: lkf_vec
	@echo "--- Running Vector LKF (M4) ---"
	$(RUNNER_VEC) ./lkf_vec $(NOISY_CSV) $(LKF_VEC_OUT)

run_ekf_vec: ekf_vec
	@echo "--- Running Vector EKF (M4) ---"
	$(RUNNER_VEC) ./ekf_vec $(NOISY_CSV) $(EKF_VEC_OUT)

# Verification

verify: $(LKF_OUT) $(EKF_OUT)
	@echo "--- Verification: M3 scalar vs M2 C++ ---"
	python3 verify.py \
		$(M2_LKF_OUT) $(LKF_OUT) \
		$(M2_EKF_OUT) $(EKF_OUT)

verify_vec: $(LKF_VEC_OUT) $(EKF_VEC_OUT)
	@echo "--- Verification: M4 vector vs M3 scalar ---"
	python3 verify.py \
		$(LKF_OUT) $(LKF_VEC_OUT) \
		$(EKF_OUT) $(EKF_VEC_OUT)

# Clean

clean:
	rm -f *.o
	rm -f lkf_asm ekf_asm lkf_vec ekf_vec
	rm -f $(LKF_OUT) $(EKF_OUT) $(LKF_VEC_OUT) $(EKF_VEC_OUT)
	@echo "  Cleaned."

help:
	@echo ""
	@echo "  Milestone 3 (scalar):"
	@echo "    make all          build lkf_asm and ekf_asm"
	@echo "    make run_lkf      run scalar LKF"
	@echo "    make run_ekf      run scalar EKF"
	@echo "    make verify       compare M3 vs M2"
	@echo ""
	@echo "  Milestone 4 (vector):"
	@echo "    make all_vec      build lkf_vec and ekf_vec"
	@echo "    make run_lkf_vec  run vector LKF"
	@echo "    make run_ekf_vec  run vector EKF"
	@echo "    make verify_vec   compare M4 vs M3"
	@echo ""
	@echo "    make all_m4       build everything"
	@echo "    make clean        remove all outputs"