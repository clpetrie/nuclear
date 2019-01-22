set terminal gif
set output 'data.gif'
#plot 'data.dat' using 5:($7) pt 7 with yerrorbars
plot 'data.dat' using 5 pt 7
