linearpt=5
ippt=7
quadpt=13
exppt=3
myps=2.5
myshift=0.6
set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'energies.eps'
set xlabel "System"
set ylabel "Energy per nucleon (MeV)"
set xrange [2:30]
set xtics ("{}^4He" 4, "{}^{16}O" 16, "SNM" 28) #offset 0,-0.3
#set key at 22, -35 box lw 2
set key at 20, -14 box lw 2
plot 'energy.dat' using ($1-1.5*myshift):2:3 with errorbars title 'Linear' lc rgb 'blue' ps myps pt linearpt,\
   'energy.dat' using ($1-0.5*myshift):4:5 with errorbars title 'Independent Pair' lc rgb 'red' ps myps pt ippt,\
   'energy.dat' using ($1+0.5*myshift):6:7 with errorbars title 'Quadratic' lc rgb 'green' ps myps pt quadpt,\
   'energy.dat' using ($1+1.5*myshift):8 title 'Experimental' lc rgb 'orange' ps myps pt exppt
