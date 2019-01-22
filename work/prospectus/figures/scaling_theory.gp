set terminal postscript eps enhanced color font 'Helvetica, 25'
set terminal postscript eps color
set output 'scaling_theory.eps'
linewidth=5
xmax=28

set xlabel "A"
set ylabel "Number of Terms"
flin(x)=x*(x-1)/2
fip(x)=x*(x-1)*(x-2)*(x-3)/8
fquad(x)=flin(x)**2-flin(x)
sip(x)=(fip(x)+flin(x))/flin(x)
squad(x)=(fquad(x)+flin(x))/flin(x)
set key at xmax-xmax*0.25, squad(xmax)-squad(xmax)*0.1 box lw 2
set xrange [1:xmax]
set yrange [0:]
set xtics (4 4, 8 8, 12 12, 16 16, 20 20, 24 24, 28 28)
#plot flin(x) title 'Linear' dashtype 1 lc rgb 'blue' lw linewidth, fip(x) title 'Independent Pair' dashtype 2 lc rgb 'red' lw linewidth+6, fquad(x) title 'Quadratic' dashtype 3 lc rgb 'green' lw linewidth+6
plot sip(x) title 'Independent Pair' dashtype 2 lc rgb 'blue' lw linewidth+6, squad(x) title 'Quadratic' dashtype 3 lc rgb 'red' lw linewidth+6
