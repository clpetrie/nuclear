set terminal postscript eps enhanced color font 'Helvetica, 25'
set output 'scaling.eps'
set xlabel "System"
set ylabel "Scaling Factor"
set xrange [2:42]
set yrange [1:1500]
set xtics ("{}^4He" 4, "{}^{16}O" 16, "{}^{28}SNM" 28, "{}^{40}Ca" 40) #offset 0,-0.3
set key at 30, 1000 box lw 2
set style func linespoints
plot 'scaling.dat' using 1:2 title 'IP Quadratic' lc rgb 'blue' ps 2.5 pt 7, \
   'scaling.dat' using 1:3 title 'Quadratic' lc rgb 'red' ps 1.5 pt 5
