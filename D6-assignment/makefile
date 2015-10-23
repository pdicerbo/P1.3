CC=gcc
CFLAGS=
LIBS= -lm

all: Heat_equation_template.x 

Heat_equation_template.x: %.x: %.c
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@ 

clean: 
	rm -f *.x
