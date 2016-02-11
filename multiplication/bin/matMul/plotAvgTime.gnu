#!/usr/bin/gnuplot

# Gnuplot script file for plotting data in file "timing.txt"
# This file is called   plotScript.gnu
# tutorials at http://www.duke.edu/~hpgavin/gnuplot.html
# http://linux.byexamples.com/archives/487/plot-your-graphs-with-command-line-gnuplot/

set terminal png
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Average Performance analysis of matrix-matrix multiplication"
set xlabel "Size of the vector."
set ylabel "Average time (sec)"
set style data linespoints
set output 'timingMatMulAvg.png'
plot  "avgMatMul.txt" using 1:2  title "Average Time"
