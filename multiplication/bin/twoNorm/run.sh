#!/bin/bash

echo Executing 2 norm of vector program!
for i in $(seq 1 10)
do
    filename=$(printf "time_%d.txt" "$i")
    for j in $(seq 10 10 1000)
    do
        ./twoNorm $j >> $filename
    done
done
echo Producing Performance Analysis!
	./avg.py
	./plotTime.gnu
	./plotFlops.gnu
	./plotAvgTime.gnu
	./plotAvgFlops.gnu
echo 2-norm of vector program Finished!

