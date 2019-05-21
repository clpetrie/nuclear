filename='potential_4he.out'
set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'vij.eps'
mylw=5.0
set xlabel "r (fm)"
set ylabel "v_n(r)"
set xrange [0.0:2.0]
set key at 1.95, 95
#set xzeroaxis
plot filename using 1:3 with lines title 'v_2' lc rgb 'blue' lw mylw,\
   filename using 1:4 with lines title 'v_3' lc rgb 'red' lw mylw,\
   filename using 1:5 with lines title 'v_4' lc rgb 'green' lw mylw,\
   filename using 1:6 with lines title 'v_5' lc rgb 'black' lw mylw,\
   filename using 1:7 with lines title 'v_6' lc rgb 'orange' lw mylw
