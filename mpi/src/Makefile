CC=mpicc
FLAG=-g -Wall
INCLUDE=-I/usr/lib/openmpi/include

all: ring grid pipe pdgemm

ring: ring_mpi.c
	${CC} ${FLAG} ${INCLUDE} ring_mpi.c -o ring

grid: grid_mpi.c
	${CC} ${FLAG} ${INCLUDE} grid_mpi.c -o grid

pipe: pipe_mpi.c
	${CC} ${FLAG} ${INCLUDE} pipe_mpi.c -o pipe

pdgemm: pdgemm.c
	${CC} ${FLAG} ${INCLUDE} pdgemm.c -o pdgemm -lcblas

clean:
	@rm -f *.o
	@rm -fr *.dSYM
	@rm -f *~

pristine: clean
	@rm -f ring
	@rm -f grid
	@rm -f pipe
	@rm -f pdgemm
