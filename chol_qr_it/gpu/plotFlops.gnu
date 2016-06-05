#!/usr/bin/gnuplot

# Gnuplot script file for plotting data in file "timing.txt"
# This file is called   plotScript.gnu
# tutorials at http://www.duke.edu/~hpgavin/gnuplot.html
# http://linux.byexamples.com/archives/487/plot-your-graphs-with-command-line-gnuplot/

set terminal png
unset label                            # remove any previous labels
set xtics ('2048' 1,'4096' 2,'8192' 3,'16384' 4,'32768' 5,'65536' 6,'131072' 7) 
set ytic auto                          # set ytics automatically
set xlabel "Size of the matrix(*128)"
set ylabel "Flops (Gflops/sec)"
set style data linespoints
set output 'flops.png'
plot  "flops.txt" using 1:2  title "lapack" ,\
      "flops.txt" using 1:3  title "CPU" ,\
      "flops.txt" using 1:4  title "GPU"
