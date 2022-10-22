# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
CCX=g++

# main executables 
EXECS = optsais optsais64 permute

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

lib/optsais32.o: lib/optSAIS.cpp lib/optSAIS.h
	$(CCX) $(CXX_FLAGS) -c -o $@ $< 

lib/optsais64.o: lib/optSAIS.cpp lib/optSAIS.h
	$(CCX) $(CXX_FLAGS) -c -o $@ $< -DM64

optsais: main.cpp common.hpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais32.o 
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp common.hpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais32.o -ldl

optsais64: main.cpp common.hpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais64.o 
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp common.hpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais64.o  -ldl -DM64 

permute: permuteBWT.cpp external/malloc_count/malloc_count.o 
	$(CXX) $(CXX_FLAGS) -o $@ permuteBWT.cpp external/malloc_count/malloc_count.o -ldl #-fsanitize=address

clean:
	rm -f $(EXECS) $(EXECS_NT) lib/*.o  