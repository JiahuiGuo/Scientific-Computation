#!/bin/bash

echo Executing omp_trsm!
for i in $(seq 1 10)
do
    filename=$(printf "time_%d.txt" "$i")
    for j in $(seq 100 100 3000)
    do
        ./omp_trsm $j $j >> $filename
    done
done
echo Producing Performance Analysis!
	./plotTime.gnu
	./avg.py
	./plotAvgTime.gnu
echo Omp_trsm program finished!

