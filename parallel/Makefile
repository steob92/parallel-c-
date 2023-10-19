# Compiler
CC=g++
NV=nvcc
#CFlags
CFLAGS=-fopenmp
NVFLAGS=-Xcompiler $(CFLAGS)

OBJECTS = calculate_pi calculate_pi_parallel 
CUOBJECTS=calculate_pi_gpu
all:$(OBJECTS) $(CUOBJECTS)

$(OBJECTS): %: %.cpp
	$(CC) $(CFLAGS)  $< -o $@ 

$(CUOBJECTS): %: %.cu
	$(NV) $(NVFLAGS)  $< -o $@ 

clean:
	rm -rf $(OBJECTS)
