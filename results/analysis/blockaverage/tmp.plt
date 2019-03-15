set terminal png enhanced font 'Verdana,33' size 1920,1440
set output 'blockaverage.png'
set xlabel 'Block Size' font "Times-Roman,54"
set ylabel 'Statistical Error' font "Times-Roman,54" offset 2.5,0
set xrange[0:85]
plot "tmp.dat" i 2 w lp lw 6 ps 4.5 notitle
