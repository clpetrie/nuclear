set term png font "Helvetica,15"
set output 'plot.png'
set xlabel '# aux. fields'
set ylabel 'E0 (MeV)'
set key off
#set title "Alpha Ground State Energy"
set style line 1 lc rgb '#dd181f' lt 1 lw 2 pt 7 ps 1.0
plot 'plotdata.txt' using 1:2:3 w yerrorbars ls 1, '' using 1:2 w lines ls 1
