linpt=5
ippt=7
myps=2.5
idmc=1 #0=vmc 1=dmc
set terminal postscript eps enhanced color font 'Helvetica, 23'
set output 'plot.eps'
set xlabel "{/Symbol r} (fm^{-3})" enhanced
set ylabel "E/A (MeV)"
set xrange [0.148:0.172]
#set yrange [-16.5:-7.0] #for vmc
set yrange [-16.5:-13.0] #for dmc
#set key at 22, -35 box lw 2
#set key at 3.0, -10 box lw 2
#set object 5 rect from 0.157 -15.29 to 0.171, -16.43 fc rgb "cyan" fs pattern 1 bo -1 #box for saturation propterties
#set object 5 rect from 0.160 -15.29 to 0.165, -16.43 fc rgb "cyan" fs pattern 1 bo -1 #box for saturation propterties
#set object 1 rect from graph 0.165, graph -9 to graph 0.160, graph -8 back
#set object 1 rect fc rgb "cyan" fillstyle solid 1.0
#set object 2 rect from 0.160 ,-9 to 0.165,-8 fc rgb 'blue' lt 1
#set object 2 rect default
set object 1 rect from 0.157 ,-16.43 to 0.171,-15.29 fc lt 3
set object 1 rect fc rgb "grey" fillstyle solid 0.5
plot 'data.txt' using 1:2+2*idmc:3+2*idmc with errorbars title 'Linear' lc rgb 'blue' ps myps pt linpt,\
   'data.txt' using 1:6+2*idmc:7+2*idmc with errorbars title 'IP Quad' lc rgb 'red' ps myps pt ippt
