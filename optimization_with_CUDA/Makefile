# Define compilers
CXX = g++
NVCC = nvcc

# Define compiler flags
CXXFLAGS = -Wall -Wextra -std=c++17 -DNDEBUG -O3
NVCCFLAGS = #-arch=sm_75 # Adjust compute capability for your GPU
LDFLAGS = -L/usr/local/cuda/lib64 -lcudart -lcuda -lcublas -lcusolver -lgtest -lgtest_main # Link CUDA runtime library

# Source files
CPP_SOURCES = geoopt_cuda_unittest.cc
CUDA_SOURCES = geoopt_cuda.cu

# Object files
CPP_OBJECTS = $(CPP_SOURCES:.cc=.o)
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=.o)

# Executable name
TARGET = geoopt_cuda_unittest

all: $(TARGET)

$(TARGET): $(CPP_OBJECTS) $(CUDA_OBJECTS)
	$(CXX) $(CPP_OBJECTS) $(CUDA_OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(CPP_OBJECTS) $(CUDA_OBJECTS) $(TARGET)
