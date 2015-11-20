#CC = gcc
CC = nvcc
LIBS = 
CFLAGS = 
EXEC = cuda_stuff.x

all: $(EXEC)

# omp: LIBS += -fopenmp
# omp: CFLAGS += -fopenmp -DLINK_OMP
# omp: all

%.x : %.o
	$(CC) $^ $(LIBS) -o $@

%.o : %.cu
	$(CC) $(CFLAGS) -c $^

clean:
	rm -rf *.o *.x