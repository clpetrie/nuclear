set terminal pngcairo size 1024,2160
set output 'plots.png'
set lmargin at screen 0.15
set rmargin at screen 0.95
set multiplot layout 3,1 scale 1,0.333
set size 1,0.33333
set origin 0,0.66666
set title 'Observable and average'
plot 'tmp.dat' i 0 u :2:3 w err t 'Observable', \
     'tmp.dat' i 0 u :4:5 w err t 'Average'
set size 1,0.33333
set origin 0,0.33333
set title 'Autocorrelations'
unset key
plot 'tmp.dat' i 1 w lp t 'Autocorrelations'
set size 1,0.33333
set origin 0,0
set title 'Errors as a function of block size'
plot 'tmp.dat' i 2 w lp t 'Errors as a function of block size'
