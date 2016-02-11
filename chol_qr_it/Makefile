CC=gcc
LIBS=-L/usr/local/lib/
INCLUDE=-I../include
CCFLAGS=-g -O3 -Wall ${INCLUDE}
OBJECT=chol_qr_it.o

all:${OBJECT}
	${CC} ${CCFLAGS} ${LIBS} ${OBJECT} -o chol_qr_it -lm /usr/local/lib/libf2c.a /usr/local/lib/libblas.a /usr/local/lib/liblapack.a /usr/local/lib/libtmglib.a -lcblas

%.o:%.[ch]
	${CC} -g -Wall ${INCLUDE} $< -c

clean:
	rm -f *.o
