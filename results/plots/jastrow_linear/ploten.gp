jastrowpt=5
ippt=7
linearpt=13
exppt=9
myps=2.5
myshift=0.1
set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'energy_jaslin.eps'
set xlabel "System"
set ylabel "Energy per nucleon (MeV)"
set xrange [0.5:3.5]
set xtics ("{}^4He" 1, "{}^{16}O" 2, "{}^{40}Ca" 3) #offset 0,-0.3
#set key at 22, -35 box lw 2
set key at 1.8, -7.80 box lw 2
#plot 'energy.dat' using ($1-1.5*myshift):3:4 with errorbars title 'Jastrow' lc rgb 'blue' ps myps pt jastrowpt,\
#   'energy.dat' using ($1-0.5*myshift):5:6 with errorbars title 'Linear' lc rgb 'red' ps myps pt linearpt,\
#   'energy.dat' using ($1-0.5*myshift):7 with errorbars title 'Exp' lc rgb 'yellow' ps myps pt exppt
plot 'energy.dat' using ($1-myshift):3:4 with errorbars title 'Jastrow' lc rgb 'blue' ps myps pt jastrowpt,\
   'energy.dat' using ($1):5:6 with errorbars title 'Linear' lc rgb 'red' ps myps pt linearpt,\
   'energy.dat' using ($1+myshift):7 title 'Exp' lc rgb 'orange' ps myps pt exppt
