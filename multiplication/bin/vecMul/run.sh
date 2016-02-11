#!/bin/bash

echo Executing vector and matrix multiplication program!
for i in $(seq 1 10)
do
	filename=$(printf "time_%d.txt" "$i")
	for j in $(seq 10 10 1000)
	do
    	./vecMul $j >> $filename
	done
done
echo Producing Performance Analysis!
	./avg.py
	./plotTime.gnu
	./plotFlops.gnu
	./plotAvgTime.gnu
	./plotAvgFlops.gnu
echo Vector and matrix multiplication program Finished!
