reset

###terminal setup
#set terminal epslatex color
#set output 'fig.tex'
set terminal x11 enhanced size 600,800

###start plots
set multiplot layout 2,1

###plot average energies
set origin 0,0.5
#set xlabel "Averaging block"
set ylabel "E (MeV)"
set key outside
plot file1 using 4 with lines title "no ind pair", file2 using 4 with lines title "ind pair"

###plot errors
set origin 0,0
set xlabel "Averaging block"
set ylabel "{/Symbol D}E (MeV)"
plot file1 using 6 with lines title "no ind pair", file2 using 6 with lines title "ind pair"

unset multiplot
#pause -1
