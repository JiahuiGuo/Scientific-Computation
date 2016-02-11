all:	
	@cd src; make all

install: all
	@cd src; make install

run: install
	@cd pdf; make pdf

report: run
	@cd bin; ./run.sh

release: report
	@make pristine
	tar cvfz JiahuiGuo.tar include src bin Makefile pdf

clean:
	@rm -f *.tar
	@cd src; make clean
	@cd pdf; make clean

pristine: clean
	@cd src; make pristine
	@cd pdf; make clean
	@cd bin; rm -f ./twoNorm/twoNorm ./vecMul/vecMul ./matMul/matMul ./twoNorm/*.txt ./vecMul/*.txt ./matMul/*.txt ./twoNorm/*.png ./vecMul/*.png ./matMul/*.png  


