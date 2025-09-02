# Compiler
CXX = g++

CXXFLAGS := -I. -Wall -Wextra -std=c++20
CXXFLAGS += -DNDEBUG -O3 -pthread -march=native -flto

.PHONY: all
MKL_ROOT := /opt/intel/oneapi/mkl/2024.2
RPATH := -Wl,-rpath,$(MKL_ROOT)/lib

LDFLAGS := -L$(MKL_ROOT)/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lgtest -lgtest_main -lpthread
CXXFLAGS += -I$(MKL_ROOT)/include

# Source files
SRCS = geoopt_unittest.cc geoopt.cc

# Object files (generated from source files)
OBJS = $(SRCS:.cc=.o)

# Executable name
TARGET = geoopt_unittest

# Default rule: build the executable
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) $(RPATH)

# Rule for compiling .cc files into .o files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Phony target for cleaning up generated files
.PHONY: clean

clean:
	rm -f $(OBJS) $(TARGET)
