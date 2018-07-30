# Compiler
CC=g++-8

# Any compiler flags
CFLAG = -std=c++11 -O2

# Build Target Executable
TARGET = LidarMC++

# LINKER
LINKER = -larmadillo

# top-level rule, to compile everything.
all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAG) $(LINKER) -o $(TARGET) $(TARGET).cpp
