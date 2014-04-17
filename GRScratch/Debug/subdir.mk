################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Function_PCG_WoodScratch.c \
../blas_malloc.c \
../blas_sparse_L1_double.c \
../blas_sparse_L23_double.c \
../blas_sparse_handle.c \
../csr_double.c \
../microblas_double.c \
../spblasi_matrix.c \
../spblasi_matrix_double.c \
../spblasi_matrix_prop.c \
../spblasi_table.c \
../spvec_double.c \
../table.c \
../testRun_Function_PCG_WoodScratch.c 

OBJS += \
./Function_PCG_WoodScratch.o \
./blas_malloc.o \
./blas_sparse_L1_double.o \
./blas_sparse_L23_double.o \
./blas_sparse_handle.o \
./csr_double.o \
./microblas_double.o \
./spblasi_matrix.o \
./spblasi_matrix_double.o \
./spblasi_matrix_prop.o \
./spblasi_table.o \
./spvec_double.o \
./table.o \
./testRun_Function_PCG_WoodScratch.o 

C_DEPS += \
./Function_PCG_WoodScratch.d \
./blas_malloc.d \
./blas_sparse_L1_double.d \
./blas_sparse_L23_double.d \
./blas_sparse_handle.d \
./csr_double.d \
./microblas_double.d \
./spblasi_matrix.d \
./spblasi_matrix_double.d \
./spblasi_matrix_prop.d \
./spblasi_table.d \
./spvec_double.d \
./table.d \
./testRun_Function_PCG_WoodScratch.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


