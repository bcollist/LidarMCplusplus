# Compiler
CC=g++-7

# Any compiler flags
CFLAG = -std=c++11 -O2

# Build Target Executable
TARGET = LidarMC++

# LINKER
LINKER = -larmadillo


#INCLUDE = $ARMADILLO_ROOT/include
#LINKER = $ARMADILLO_ROOT/LIB -larmadillo


# top-level rule, to compile everything.
all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAG) $(LINKER) -o $(TARGET) $(TARGET).cpp
