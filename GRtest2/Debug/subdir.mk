################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Function_PCG_Wood_clean_codegen2.c \
../Function_PCG_Wood_clean_codegen2_emxutil.c \
../Function_PCG_Wood_clean_codegen2_initialize.c \
../Function_PCG_Wood_clean_codegen2_terminate.c \
../abs.c \
../colon.c \
../cond.c \
../diag.c \
../eye.c \
../mrdivide.c \
../rdivide.c \
../rtGetInf.c \
../rtGetNaN.c \
../rt_nonfinite.c \
../sum.c \
../testRun_Function_PCG_Wood_clean_codegen.c 

OBJS += \
./Function_PCG_Wood_clean_codegen2.o \
./Function_PCG_Wood_clean_codegen2_emxutil.o \
./Function_PCG_Wood_clean_codegen2_initialize.o \
./Function_PCG_Wood_clean_codegen2_terminate.o \
./abs.o \
./colon.o \
./cond.o \
./diag.o \
./eye.o \
./mrdivide.o \
./rdivide.o \
./rtGetInf.o \
./rtGetNaN.o \
./rt_nonfinite.o \
./sum.o \
./testRun_Function_PCG_Wood_clean_codegen.o 

C_DEPS += \
./Function_PCG_Wood_clean_codegen2.d \
./Function_PCG_Wood_clean_codegen2_emxutil.d \
./Function_PCG_Wood_clean_codegen2_initialize.d \
./Function_PCG_Wood_clean_codegen2_terminate.d \
./abs.d \
./colon.d \
./cond.d \
./diag.d \
./eye.d \
./mrdivide.d \
./rdivide.d \
./rtGetInf.d \
./rtGetNaN.d \
./rt_nonfinite.d \
./sum.d \
./testRun_Function_PCG_Wood_clean_codegen.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


