#!/bin/bash

# run 2-Norm of vector
cd ./twoNorm
./run.sh
cd ..
# copy the png files for report
cp ./twoNorm/*.png ../pdf/

#run vector matrix multiplication
cd ./vecMul
./run.sh
cd ..
# copy the png files for report
cp ./vecMul/*.png ../pdf/

#run matrix matrix multiplication
cd ./matMul
./run.sh
cd ..
# copy the png files for report
cp ./matMul/*.png ../pdf/

#run to generate the average plot
./plotAvgTime.gnu
./plotAvgFlops.gnu
cp *.png ../pdf/
