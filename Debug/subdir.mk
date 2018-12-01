################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../MMFF94Parameters.cpp \
../forcefieldmmff94.cpp \
../main.cpp 

OBJS += \
./MMFF94Parameters.o \
./forcefieldmmff94.o \
./main.o 

CPP_DEPS += \
./MMFF94Parameters.d \
./forcefieldmmff94.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include/openbabel-2.0 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


