# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -std=c++17
LDFLAGS = -lgtest -lgtest_main -lpthread

# Source files
SRCS = codewars_unittest.cc \
       bouncing_balls.cc \
       consecutive_strings.cc \
       convert_string_to_camel_case.cc \
       counting_duplicates.cc \
       detect_pangram.cc \
       duplicate_encoder.cc \
       every_possible_sum_of_two_digits.cc \
       multiples_of_3_or_5.cc

# Object files (generated from source files)
OBJS = $(SRCS:.cc=.o)

# Executable name
TARGET = codewars_unittest

# Default rule: build the executable
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Rule for compiling .cc files into .o files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Phony target for cleaning up generated files
.PHONY: clean

clean:
	rm -f $(OBJS) $(TARGET)
