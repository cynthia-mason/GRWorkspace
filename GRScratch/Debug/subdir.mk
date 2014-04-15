################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Function_PCG_WoodScratch.c \
../testRun_Function_PCG_WoodScratch.c 

OBJS += \
./Function_PCG_WoodScratch.o \
./testRun_Function_PCG_WoodScratch.o 

C_DEPS += \
./Function_PCG_WoodScratch.d \
./testRun_Function_PCG_WoodScratch.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


