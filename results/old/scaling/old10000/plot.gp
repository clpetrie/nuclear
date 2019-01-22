set terminal postscript eps size 3.5,2.62 enhanced color \
    font 'Helvetica,20' linewidth 2
set output 'scaling.eps'
set xtics 4000/4
set xlabel "# of cores"
set ylabel "1/t (sec^{-1})"
slope=0.00013958204946927416
plot 'prop.out' using ($1):(1/$2) with points pointtype 7 lc rgb "blue" notitle, \
   slope*x lc rgb "black" notitle
