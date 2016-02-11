#!/usr/bin/gnuplot

# Gnuplot script file for plotting data in file "timing.txt"
# This file is called   plotScript.gnu
# tutorials at http://www.duke.edu/~hpgavin/gnuplot.html
# http://linux.byexamples.com/archives/487/plot-your-graphs-with-command-line-gnuplot/

set terminal png
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Comparison of execution time of different chunk size"
set xlabel "Size of matrix."
set ylabel "Time (sec)"
set style data linespoints
set key left top
set output 'static_8.png'
plot  "static_ck50_8.txt" using 1:2  title "Chunk size 50" ,\
	  "static_ck10_8.txt" using 1:2  title "Chunk size 10"
