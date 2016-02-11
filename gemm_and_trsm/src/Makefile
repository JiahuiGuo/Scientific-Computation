VPATH=./src

CC=gcc
CCFLAG=-g -Wall -fopenmp

all: omp_gemm omp_trsm
	
omp_gemm:omp_gemm.c
	${CC} ${CCFLAG} omp_gemm.c c_timer.c -o omp_gemm

omp_trsm:omp_trsm.c
	${CC} ${CCFLAG} omp_trsm.c c_timer.c -o omp_trsm

install:all
	@cp -r omp_gemm ../bin/omp_gemm/
	@cp -r omp_trsm ../bin/omp_trsm/

clean:
	@rm -f *.o

pristine: clean
	@rm -f omp_gemm omp_trsm
