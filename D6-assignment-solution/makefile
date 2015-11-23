CC=mpicc
CFLAGS=
LIBS= -lm

all: Heat_equation_mpi.x

Heat_equation_mpi.x: %.x: %.c
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

clean: 
	rm -f *.x
