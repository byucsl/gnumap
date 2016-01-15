
OBJS=align.cu.o align.cpp.o
EXE=driver

#FLAGS=-deviceemu -g -DDEBUG
#FLAGS=-deviceemu -g
#FLAGS=-DSEQ_STYLE
FLAGS=-O2

all: $(OBJS) driver.cpp
	nvcc $(FLAGS) -o $(EXE) driver.cpp $(OBJS)


align.cpp.o: align.cpp
	nvcc $(FLAGS) -o align.cpp.o -c align.cpp

align.cu.o: align.cu
	nvcc $(FLAGS) -o align.cu.o -c align.cu

clean:
	rm -f *.o
	rm -f $(EXE)
