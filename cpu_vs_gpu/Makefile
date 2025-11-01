# Define compilers
CXX = g++
NVCC = nvcc

# Define compiler flags
CXXFLAGS = -Wall -Wextra -std=c++17 -DNDEBUG -O3 -fopenmp
NVCCFLAGS = #-arch=sm_75 # Adjust compute capability for your GPU
LDFLAGS = -L/usr/local/cuda/lib64 -lcudart -lcuda -lgtest -lgtest_main -fopenmp # Link CUDA runtime library

# Source files
CPP_SOURCES = integration_cpu.cc cpu_vs_gpu_unittest.cc
CUDA_SOURCES = integration_gpu.cu

# Object files
CPP_OBJECTS = $(CPP_SOURCES:.cc=.o)
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=.o)

# Executable name
TARGET = cpu_vs_gpu_unittest

all: $(TARGET)

$(TARGET): $(CPP_OBJECTS) $(CUDA_OBJECTS)
	$(CXX) $(CPP_OBJECTS) $(CUDA_OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(CPP_OBJECTS) $(CUDA_OBJECTS) $(TARGET)
