linearpt=5
ippt=7
quadpt=13
exppt=9
myps=2.5
myshift=0.1
set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'energy.eps'
set xlabel "System"
set ylabel "Energy per nucleon (MeV)"
set xrange [0.5:4.5]
#set xtics ("{}^4He" 4, "{}^{16}O" 16, "SNM" 28) #offset 0,-0.3
set xtics ("{}^4He" 1, "{}^{16}O" 2, "{}^{40}Ca" 3, "SNM" 4) #offset 0,-0.3
#set key at 22, -35 box lw 2
set key at 3.0, -10 box lw 2
plot 'energy.dat' using ($1-1.5*myshift):3:4 with errorbars title 'Linear' lc rgb 'blue' ps myps pt linearpt,\
   'energy.dat' using ($1-0.5*myshift):5:6 with errorbars title 'IP Quadratic' lc rgb 'red' ps myps pt ippt,\
   'energy.dat' using ($1+0.5*myshift):7:8 with errorbars title 'Quadratic' lc rgb 'green' ps myps pt quadpt,\
   'energy.dat' using ($1+1.5*myshift):9 title 'Experimental' lc rgb 'orange' ps myps pt exppt
