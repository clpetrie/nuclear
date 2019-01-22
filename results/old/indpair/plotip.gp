linearpt=5
ippt=7
exppt=3
myps=2.5
myshift=0.6
set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'energiesip.eps'
set xlabel "System"
set ylabel "Energy per nucleon (MeV)"
set xrange [2:30]
set xtics ("{}^4He" 4, "{}^{16}O" 16, "SNM" 28) #offset 0,-0.3
#set key at 22, -35 box lw 2
set key at 20, -14 box lw 2
plot 'energy.dat' using ($1-1*myshift):2:3 with errorbars title 'Linear' lc rgb 'blue' ps myps pt linearpt,\
   'energy.dat' using ($1-0*myshift):4:5 with errorbars title 'Independent Pair' lc rgb 'red' ps myps pt ippt,\
   'energy.dat' using ($1+1*myshift):8 title 'Experimental' lc rgb 'orange' ps myps pt exppt
