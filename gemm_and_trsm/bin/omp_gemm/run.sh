#!/bin/bash

echo Executing omp_gemm!
for i in $(seq 1 10)
do
    filename=$(printf "time_%d.txt" "$i")
    for j in $(seq 100 100 3000)
    do
        ./omp_gemm $j $j $j >> $filename
    done
done
echo Producing Performance Analysis!
	./plotTime.gnu
	./avg.py
echo Omp_gemm program finished!

