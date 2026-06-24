

# Toolchain
CC      = riscv64-linux-gnu-gcc
AS      = riscv64-linux-gnu-gcc    
LD      = riscv64-linux-gnu-gcc
RUNNER  = qemu-riscv64                 

# Compiler flags
CFLAGS  = -march=rv64imafd -mabi=lp64d -O0 -g -Wall
ASFLAGS = -march=rv64imafd -mabi=lp64d -g
LDFLAGS = -march=rv64imafd -mabi=lp64d -static -lm

# Source files
ASM_COMMON = kalman_common.s
ASM_LKF    = lkf_asm.s
ASM_EKF    = ekf_asm.s
DRV_LKF    = driver_lkf.c
DRV_EKF    = driver_ekf.c

# Object files
OBJ_COMMON = kalman_common.o
OBJ_LKF    = lkf_asm.o driver_lkf.o
OBJ_EKF    = ekf_asm.o lkf_asm.o driver_ekf.o

# Data files (adjust path if needed)
NOISY_CSV  = noisy_data.csv
M2_LKF_OUT = lkf_output.csv        
M2_EKF_OUT = ekf_output_final.csv   

# Output files 
LKF_OUT    = lkf_asm_output.csv
EKF_OUT    = ekf_asm_output.csv

# Build targets

.PHONY: all clean run_lkf run_ekf verify help

all: lkf_asm ekf_asm
	@echo ""
	@echo "  Build complete."
	@echo "  Run:  make run_lkf    to execute LKF"
	@echo "  Run:  make run_ekf    to execute EKF"
	@echo "  Run:  make verify     to check against Milestone 2 output"

# Compile common assembly
$(OBJ_COMMON): $(ASM_COMMON)
	$(AS) $(ASFLAGS) -c $< -o $@

# Compile LKF assembly
lkf_asm.o: $(ASM_LKF)
	$(AS) $(ASFLAGS) -c $< -o $@

# Compile LKF C driver
driver_lkf.o: $(DRV_LKF)
	$(CC) $(CFLAGS) -c $< -o $@

# Link LKF executable
lkf_asm: $(OBJ_COMMON) lkf_asm.o driver_lkf.o
	$(LD) $(LDFLAGS) $^ -o $@
	@echo "  Linked: lkf_asm"

# Compile EKF assembly
ekf_asm.o: $(ASM_EKF)
	$(AS) $(ASFLAGS) -c $< -o $@

# Compile EKF C driver 
driver_ekf.o: $(DRV_EKF)
	$(CC) $(CFLAGS) -c $< -o $@

# Link EKF executable 
ekf_asm: $(OBJ_COMMON) $(OBJ_EKF)
	$(LD) $(LDFLAGS) $^ -o $@
	@echo "  Linked: ekf_asm"

# Run targets

run_lkf: lkf_asm
	@echo "--- Running LKF Assembly ---"
	$(RUNNER) ./lkf_asm $(NOISY_CSV) $(LKF_OUT)

run_ekf: ekf_asm
	@echo "--- Running EKF Assembly ---"
	$(RUNNER) ./ekf_asm $(NOISY_CSV) $(EKF_OUT)

# Verification

verify: $(LKF_OUT) $(EKF_OUT)
	@echo "--- Numerical Verification (M3 vs M2) ---"
	python3 verify.py \
		$(M2_LKF_OUT) $(LKF_OUT) \
		$(M2_EKF_OUT) $(EKF_OUT)


# Clean
clean:
	rm -f *.o lkf_asm ekf_asm $(LKF_OUT) $(EKF_OUT)
	@echo "  Cleaned."

help:
	@echo "Targets: all  run_lkf  run_ekf  verify  clean"