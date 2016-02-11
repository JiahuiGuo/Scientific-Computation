VPATH=./src

TARGET_TwoNorm=twoNorm
OBJECTS_TwoNorm=twoNorm.o c_timer.o

TARGET_VecMul=vecMul
OBJECTS_VecMul=vecMul.o c_timer.o

TARGET_MatMul=matMul
OBJECTS_MatMul=matMul.o c_timer.o

CC=gcc
CCFLAG=-g -Wall

${TARGET_TwoNorm}:${OBJECTS_TwoNorm}
	${CC} ${CCFLAG} ${OBJECTS_TwoNorm} -o ${TARGET_TwoNorm} -lm -lcblas
	
${TARGET_VecMul}:${OBJECTS_VecMul}
	${CC} ${CCFLAG} ${OBJECTS_VecMul} -o ${TARGET_VecMul} -lm -lcblas

${TARGET_MatMul}:${OBJECTS_MatMul}
	${CC} ${CCFLAG} ${OBJECTS_MatMul} -o ${TARGET_MatMul} -lm -lcblas

%.o:%.[ch]  
	${CC} ${CCFLAG} $< -c

all: ${TARGET_TwoNorm} ${TARGET_VecMul} ${TARGET_MatMul}

install:
	@cp ${TARGET_TwoNorm} ../bin/twoNorm
	@cp ${TARGET_VecMul} ../bin/vecMul
	@cp ${TARGET_MatMul} ../bin/matMul
	
clean:
	@rm -f ${OBJECTS_TwoNorm} ${OBJECTS_VecMul} ${OBJECTS_MatMul}

pristine: clean
	@rm -f ${TARGET_TwoNorm} ${TARGET_VecMul} ${TARGET_MatMul}
