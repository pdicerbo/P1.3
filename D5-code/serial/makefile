CC=gcc
CFLAGS=
LIBS= -lm

all: Heat_equation_serial.x 

Heat_equation_serial.x: %.x: %.c
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@ 

clean: 
	rm -f *.x
