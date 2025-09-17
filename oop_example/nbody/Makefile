# Define compilers
CXX = g++

# Define compiler flags
CXXFLAGS = -Wall -Wpedantic -Wextra -Werror -Wshadow -std=c++20 -DNDEBUG -O3
LDFLAGS = -lgtest -lgtest_main

# Source files
CPP_SOURCES = nbody_unittest.cc nbody.cc

# Object files
CPP_OBJECTS = $(CPP_SOURCES:.cc=.o)

# Executable name
TARGET = nbody_unittest

all: $(TARGET)

$(TARGET): $(CPP_OBJECTS)
	$(CXX) $(CPP_OBJECTS) -o $@ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(CPP_OBJECTS) $(TARGET)
