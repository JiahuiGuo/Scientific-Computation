all:	
	@cd src; make all

install: all
	@cd src; make install

run: install
	@cd bin; ./run.sh

release: run
	@make pristine
	tar cvfz JiahuiGuo.tar include src bin Makefile hw2.pdf

clean:
	@rm -f *.tar
	@cd src; make clean

pristine: clean
	@cd src; make pristine
	@cd bin; rm -f ./omp_gemm/omp_gemm ./omp_gemm/*.txt ./omp_gemm/*.png
	@cd bin; rm -f ./omp_trsm/omp_trsm ./omp_trsm/*.txt ./omp_trsm/*.png


