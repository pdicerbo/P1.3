CC = gcc
LIBS = -lm
CFLAGS = -Wall -O3
SRCS = hello_world.c
OBJS = ${SRCS:.c=.o}
EXEC = hello_world.x Heat_equation_serial.x OptMatTransp.x Heat_equation_parallel.x

all: $(EXEC)

omp: LIBS += -fopenmp
omp: CFLAGS += -fopenmp -DLINK_OMP
omp: all

Heat_equation_parallel.x : Heat_equation_parallel.c
	mpicc $^ -o $@

%.x : %.o
	$(CC) $^ $(LIBS) -o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm -rf *.o *.x