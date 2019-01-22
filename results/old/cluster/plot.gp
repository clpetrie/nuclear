set term png font "Helvetica,15"
set output 'plot.png'
set xlabel 'r (fm)'
set ylabel 'g(r)'
set key bottom right
plot 'gnp.dmc.block' using 1:2:3 with yerrorbars lt rgb 'red' title 'np0', 'gnp.dmc.block' using 1:4:5 with yerrorbars lt rgb 'green' title 'np1', 'gnp.dmc.block' using 1:6:7 with yerrorbars lt rgb 'blue' title 'pp', 'gnp.dmc.block' using 1:8:9 with yerrorbars lt rgb 'violet' title 'nn', \
   'gnp.dmc.block' using 1:2 with lines lt rgb 'red' title '', 'gnp.dmc.block' using 1:4 with lines lt rgb 'green' title '', 'gnp.dmc.block' using 1:6 with lines lt rgb 'blue' title '', 'gnp.dmc.block' using 1:8 with lines lt rgb 'violet' title ''
