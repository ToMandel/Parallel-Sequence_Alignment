build:
	mpicxx -fopenmp -c prog.c -o prog.o
	nvcc -gencode arch=compute_75,code=sm_75 -c cudaFunc.cu -o cudaFunc.o
	mpicxx -fopenmp -o mpiCudaOpemMP  prog.o cudaFunc.o  /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./mpiCudaOpemMP

run:
	mpiexec -np 3 ./mpiCudaOpemMP <input.txt> output.txt