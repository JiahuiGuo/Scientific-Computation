#!/usr/bin/gnuplot

# Gnuplot script file for plotting data in file "timing.txt"
# This file is called   plotScript.gnu
# tutorials at http://www.duke.edu/~hpgavin/gnuplot.html
# http://linux.byexamples.com/archives/487/plot-your-graphs-with-command-line-gnuplot/

set terminal png
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Average Performance comparison"
set xlabel "Size of the vector."
set ylabel "Average Flops(Gflops/sec)"
set style data linespoints
set output 'flopsAvg.png'
plot  "./twoNorm/avgTwoNorm.txt" using 1:3  title "2-norm of vector" ,\
	  "./vecMul/avgVecMul.txt" using 1:3  title "vector * matrix" ,\
	  "./matMul/avgMatMul.txt" using 1:3  title "matrix * matrix"
