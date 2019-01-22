set terminal gif
set output 'data.gif'
set logscale y
#plot 'data.dat' using 5:($7) pt 7 with yerrorbars
plot 'data.dat' using (-$5) pt 7
#plot 'data.dat' using 5 pt 7
