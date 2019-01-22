set terminal gif enhanced
set output 'data.gif'
set title 'Correlation function calculation'
set xlabel 'm-n'
set ylabel '<phin phim>'
set style line 1 lc rgb '#8A2BE2' lt 1 lw 1.7 pt 7 ps 1.2   # --- purple
set style line 2 lc rgb '#32CD32' lt 1 lw 1.7 pt 5 ps 1.2   # --- red
plot 'data.out' using 1:2 with linespoints ls 1 title 'calculation', 'data.out' using 1:3 with linespoints ls 2 title 'theory'
