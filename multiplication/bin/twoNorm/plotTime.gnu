#!/usr/bin/gnuplot

# Gnuplot script file for plotting data in file "timing.txt"
# This file is called   plotScript.gnu
# tutorials at http://www.duke.edu/~hpgavin/gnuplot.html
# http://linux.byexamples.com/archives/487/plot-your-graphs-with-command-line-gnuplot/

set terminal png
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Performance analysis of 2 norm of vector in 10 trials"
set xlabel "Size of the vector."
set ylabel "Time (sec)"
set style data linespoints
set output 'timingTwoNorm.png'
plot  "time_1.txt" using 1:2  title "trial 1" ,\
      "time_2.txt" using 1:2  title "trial 2" ,\
      "time_3.txt" using 1:2  title "trial 3" ,\
      "time_4.txt" using 1:2  title "trial 4" ,\
      "time_5.txt" using 1:2  title "trial 5" ,\
      "time_6.txt" using 1:2  title "trial 6" ,\
      "time_7.txt" using 1:2  title "trial 7" ,\
      "time_8.txt" using 1:2  title "trial 8" ,\
      "time_9.txt" using 1:2  title "trial 9" ,\
      "time_10.txt" using 1:2  title "trial 10"
