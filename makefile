# Define output .exe file
LINK_TARGET = mc.exe

# Object Files
OBJS = \
	LidarMC++V1.o

# Compiler
CCX=clang

# Any compiler flags optimization etc.
CFLAG = -O2

INCLUDE = $ARMADILLO_ROOT/include
LINKER = $ARMADILLO_ROOT/LIB -larmadillo


# top-level rule, to compile everything.
all: $(PROG)

PROG: LidarMC++V1.o
	clang -o PROG LidarMC++V1.o

LidarMC++V1.o: LidarMC++V1.cpp
	clang -c LidarMC++V1.cpp

	clean rm LidarMC++V1.o PROG
