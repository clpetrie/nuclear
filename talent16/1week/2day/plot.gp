set terminal gif
set output 'energy.gif'
plot 'energy.dat' using 5:($7) pt 7 with yerrorbars
#plot 'energy.dat' using 5 pt 7
