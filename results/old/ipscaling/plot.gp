set terminal postscript eps enhanced color font 'Helvetica, 25'
set output 'scaling.eps'
set xlabel "System"
set ylabel "Scaling Factor"
set xrange [2:30]
set yrange [1:250]
set xtics ("{}^4He" 4, "{}^{16}O" 16, "SNM" 28) #offset 0,-0.3
#set key at 16, 230 box lw 2
plot 'scaling.dat' using 1:2 title "" lc rgb 'blue' ps 2.5 pt 7
#   'scaling.dat' using 1:3 title 'Quadratic' lc rgb 'red' ps 1.5 pt 3
