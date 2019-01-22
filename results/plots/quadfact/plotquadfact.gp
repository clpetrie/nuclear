set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'axes.tex'
set xlabel "Optimization factor $\\beta$"
set ylabel "Energy (MeV)"
set xrange [-0.01:1.01]
set key at 20, -12 box lw 2
plot 'quadfact.dat' using 1:2:3 with errorbars pt 5 ps 2 lc rgb 'blue'
