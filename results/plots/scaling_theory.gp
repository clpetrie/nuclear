set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'scaling_theory.eps'
set xlabel "Number of Particles"
set ylabel "Number of Quadratic Terms"

nip(x)=(x*(x-1)*(x-2)*(x-3))/8.0
nquad(x)=x*(x-1)/2.0*(x*(x-1)/2.0-1)
nquadfix(x)=nquad(x)-nip(x)

myps=2
a=4
b=40
set xrange [a:b]
set sample b-a+1
set key at 32.0, 500000 box lw 2
plot nip(x) w p ps myps pt 7 lc rgb 'blue' title 'IP Quadratic',\
   nquad(x) w p ps myps pt 5 lc rgb 'red' title 'Quadratic',\
   nquadfix(x) w p ps myps pt 9 lc rgb 'green' title 'Quadratic - no IP sym

#   nip with errorbars title 'IP Quadratic' lc rgb 'red'
#   'energy.dat' using ($1+0.5*myshift):7:8 with errorbars title 'Quadratic' lc rgb 'green' ps myps pt quadpt,\
#   'energy.dat' using ($1+1.5*myshift):9 title 'Experimental' lc rgb 'orange' ps myps pt exppt
