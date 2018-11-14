# Compiler
CC=g++

# Any compiler flags
CFLAG = -std=c++11

# Build Target Executable
TARGET = LidarMC++

# LINKER
LINKER = -larmadillo

# SOURCES
SOURCES = bhmie.cpp photon_tracking.cpp IOPs.cpp

# top-level rule, to compile everything.
all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAG) $(LINKER) -o $(TARGET) $(TARGET).cpp $(SOURCES)
