set terminal postscript eps enhanced color font 'Helvetica, 25'
set output 'scaling.eps'
set xlabel "# of cores"
set ylabel "1/t()"
#set xrange [2:30]
#set yrange [1:140]
set key at 24, 130 box lw 2
set style func linespoints
plot 'scaling.dat' using 1:2 title '' lc rgb 'blue' ps 2.5 pt 7
