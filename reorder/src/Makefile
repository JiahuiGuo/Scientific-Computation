VPATH=./include

CC=gcc
CCFLAG=-g -Wall
INCLUDE=-I../include
LIB=-L../lib

all: ccs reorder ccs_reorder bcrs
	
ccs:ccs.c
	${CC} ${CCFLAG} ${INCLUDE} ${LIB} ccs.c -o ccs -lpapi -lm

reorder:reorder.c
	${CC} ${CCFLAG} reorder.c -o reorder

ccs_reorder:ccs_reorder.c
	${CC} ${CCFLAG} ${INCLUDE} ${LIB} ccs_reorder.c -o ccs_reorder -lpapi -lm

bcrs:bcrs.c
	${CC} ${CCFLAG} ${INCLUDE} ${LIB} bcrs.c -o bcrs -lpapi -lm

clean:
	@rm -f *.o
	@rm -f matrix.reorder.ccs
	@rm -f matrix.reorder.crs

pristine: clean
	@rm -f ccs
	@rm -f bcrs
	@rm -f reorder
	@rm -f ccs_reorder
