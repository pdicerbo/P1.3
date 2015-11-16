CC = gcc
LIBS = -lm
CFLAGS = -Wall
SRCS = hello_world.c
OBJS = ${SRCS:.c=.o}
EXEC = hello_world.x Heat_equation_serial.x

all: $(EXEC)

omp: LIBS += -fopenmp
omp: CFLAGS += -fopenmp -DLINK_OMP
omp: all

%.x : %.o
	$(CC) $^ $(LIBS) -o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $^

clean:
	rm -rf *.o *.x