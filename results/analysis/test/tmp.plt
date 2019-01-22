set terminal x11 0 title "observable and average"
plot "tmp.dat" i 0 u :2:3 w err
replot "tmp.dat" i 0 u :4:5 w err
set terminal x11 1 title "autocorrelations"
plot "tmp.dat" i 1 w lp
set terminal x11 2 title "errors as a funcion of block size"
plot "tmp.dat" i 2 w lp
pause -1
q
