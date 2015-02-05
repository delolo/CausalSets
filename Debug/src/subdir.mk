################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CausalSets.cpp 

O_SRCS += \
../src/2d100runs.o \
../src/2d500runs.o \
../src/4d.o \
../src/4dk=40.o \
../src/MaxAndNToMax.o 

OBJS += \
./src/CausalSets.o 

CPP_DEPS += \
./src/CausalSets.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


