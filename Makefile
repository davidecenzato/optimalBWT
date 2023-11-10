# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -g
CXX=g++

# main executables 
EXECS = optsais optsais64 bcr permute 

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

lib/optsais32.o: lib/optSAIS.cpp lib/optSAIS.h
	$(CXX) $(CXX_FLAGS) -c -o $@ $< 

lib/optsais64.o: lib/optSAIS.cpp lib/optSAIS.h
	$(CXX) $(CXX_FLAGS) -c -o $@ $< -DM64

optsais: main.cpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais32.o 
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais32.o -ldl

optsais64: main.cpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais64.o 
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp computeTransform.cpp external/malloc_count/malloc_count.o lib/optsais64.o  -ldl -DM64 

bcr:
	make -C external/BCR_LCP_GSA/ SAP=1
	
permute: permuteBWT.cpp external/malloc_count/malloc_count.o 
	$(CXX) $(CXX_FLAGS) -o $@ permuteBWT.cpp external/malloc_count/malloc_count.o -ldl 

clean:
	rm -f $(EXECS) $(EXECS_NT) lib/*.o external/malloc_count/*.o
	make clean -C external/BCR_LCP_GSA/

