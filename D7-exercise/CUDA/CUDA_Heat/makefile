CC = gcc
NVCC = nvcc
LIBS = -lm
CFLAGS = 
SRC = Heat_equation_serial.cu parallel_function.cu
OBJS = ${SRC:.cu=.o}
EXEC = Heat_equation_serial.x

# omp: LIBS += -fopenmp
# omp: CFLAGS += -fopenmp -DLINK_OMP
# omp: all

cuda: CFLAGS+=-D__CUDA
cuda: parallel_function.o
cuda: OBJS+=parallel_function.o
cuda: CC=nvcc
cuda: all

all: $(EXEC)

%.x : $(OBJS)
	$(CC) $^ $(LIBS) -o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $^

%.o : %.cu
	$(CC) $(CFLAGS) -c $^

clean:
	rm -rf *.o *.x
