CCC = nvcc
CCLINKER = nvcc
CUDADIR = /data/apps/cuda/4.2/


chol_qr_it: chol_qr_main.o chol_qr_it.o
	$(CCLINKER) chol_qr_main.o  chol_qr_it.o -lm -o chol_qr_it -L$(MKLROOT)/lib/intel64 -L$(CUDADIR)/lib64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lcublas -lcudart -lstdc++ -lm -liomp5 -lgfortran  -L$(CUDADIR)/SDK/C/lib -lcutil_x86_64

chol_qr_main.o: chol_qr_main.cu
	$(CCC) chol_qr_main.cu -arch sm_20  -c -g -I$(CUDADIR)/toolkit/include -I$(CUDADIR)/SDK/C/common/inc -I/data/apps/intel/2011.5.220/composerxe-2011.5.220/mkl/include
chol_qr_it.o: chol_qr_it.cu
	$(CCC) chol_qr_it.cu -c -g -I$(CUDADIR)/toolkit/include -I$(CUDADIR)/SDK/C/common/inc -I/data/apps/intel/2011.5.220/composerxe-2011.5.220/mkl/include


clean:
	rm -f *.o
pristine: clean
	rm -f chol_qr_it
release:pristine


