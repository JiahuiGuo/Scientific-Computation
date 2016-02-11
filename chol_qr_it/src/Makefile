CC=gcc
INCLUDE=-I../include
CCFLAGS=-g -O3 -Wall ${INCLUDE}
OBJECT=chol_qr_it.o

all:${OBJECT}
	${CC} ${CCFLAGS} ${OBJECT} -o chol_qr_it -lm -lcblas ../lib/libblas.a ../lib/libf2c.a ../lib/liblapack.a ../lib/libtmglib.a

%.o:%.[ch]
	${CC} -g -Wall ${INCLUDE} $< -c

clean:
	rm -f *.o
