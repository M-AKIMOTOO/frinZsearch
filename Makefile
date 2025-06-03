# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
LDFLAGS = -lfftw3f -lfftw3f_threads -lm -lpthread -lstdc++fs

# Source files
SRCS = \
    src/frinZsearch.cpp   \
    src/frinZargs.cpp     \
    src/frinZread.cpp     \
    src/frinZlogger.cpp   \
    src/frinZparamcal.cpp \
    src/frinZfftshift.cpp \
    src/frinZfitting.cpp

# Object files (derived from SRCS)
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = frinZsearch

# Default target
all: $(TARGET)

# Rule to link the executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Rule to compile .cpp files to .o files
# This is a pattern rule that applies to all .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Target to clean up object files and the executable
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets (targets that are not actual files)
.PHONY: all clean
