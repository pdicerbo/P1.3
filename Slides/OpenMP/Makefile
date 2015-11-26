CC=icc
CFLAGS=-O3 -xAVX
CPPFLAGS=-D__OPENMP
LDFLAGS=-openmp

all: serial openmp

serial: GS_heat_openmp.c heat_openmp.c fast_transpose_openmp.c
	${CC} ${CFLAGS} GS_heat_openmp.c -o GS_heat_openmp_serial.x
	${CC} ${CFLAGS} heat_openmp.c -o heat_openmp_serial.x
	${CC} ${CFLAGS} fast_transpose_openmp.c -o fast_transpose_openmp_serial.x

openmp: GS_heat_openmp.c heat_openmp.c fast_transpose_openmp.c
	${CC} ${CFLAGS} ${CPPFLAGS} GS_heat_openmp.c -o GS_heat_openmp.x ${LDFLAGS}
	${CC} ${CFLAGS} ${CPPFLAGS} heat_openmp.c -o heat_openmp.x ${LDFLAGS}
	${CC} ${CFLAGS} ${CPPFLAGS} fast_transpose_openmp.c -o fast_transpose_openmp.x ${LDFLAGS}

clean:
	rm -rf *.x *~



