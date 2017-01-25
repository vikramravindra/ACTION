CXX=g++
INCLUDE=-I./include
CXXFLAGS=-g3 -march=native -O4 -std=c++11 -fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations -msse2  -Wall -fPIC $(INCLUDE)
SRC=src/binom.c  src/main.c src/cmp.c  src/dmvhyper.c  src/dmvhyperLog.c  src/pmvhyper.c src/CoDO.c src/pmvhyperLog.c
OBJ=$(SRC:.c=.o)
PROGRAM=CoDO_test
MEX_EXT=mexa64
	
all : $(PROGRAM) message


src/%.o: src/%.c
	$(CXX) $(CXXFLAGS) -c -o $@  $<


$(PROGRAM): $(OBJ) 
	$(CXX) $(CXXFLAGS)  -o $@ $(OBJ)	


%.mexa64: matlab/%.cpp
	mex $< $(OBJ) $(INCLUDE)

matlab: $(PROGRAM) mex_dmvhyper.mexa64 mex_pmvhyper.mexa64 mex_CoDO.mexa64


clean:
	rm -f $(PROGRAM) *.mex*  src/*.o *~


message:
	echo "Executable: $(PROGRAM) has been created"	
