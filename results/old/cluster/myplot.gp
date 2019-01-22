files='14n2p.lin_0.16 14n2p.ip_0.16'
files='14n2p.lin_0.016 14n2p.lin_0.16'
clrs='red green blue violet'
state=4
dat=2 #state=1
err=3
ttl='np0'
if (state==2) {
   dat=4;
   err=5;
   ttl='np1'
}
if (state==3) {
   dat=6;
   err=7;
   ttl='pp'
}
if (state==4) {
   dat=8;
   err=9;
   ttl='nn'
}

set term png font "Helvetica,15"
set output 'myplot.png'
set xlabel 'r (fm)'
set ylabel 'pair distribution function, g(r)'
set key bottom right
set title ttl
plot for [i=1:words(files)] word(files, i).'/gofrnp.dmc' every ::::100 using 1:dat:err with yerrorbars lt rgb word(clrs, i) title word(files, i)
#plot file1 using 1:2:3 with yerrorbars lt rgb 'red' title 'np0', 'gnp.dmc.block' using 1:4:5 with yerrorbars lt rgb 'green' title 'np1', 'gnp.dmc.block' using 1:6:7 with yerrorbars lt rgb 'blue' title 'pp', 'gnp.dmc.block' using 1:8:9 with yerrorbars lt rgb 'violet' title 'nn', \
#   'gnp.dmc.block' using 1:2 with lines lt rgb 'red' title '', 'gnp.dmc.block' using 1:4 with lines lt rgb 'green' title '', 'gnp.dmc.block' using 1:6 with lines lt rgb 'blue' title '', 'gnp.dmc.block' using 1:8 with lines lt rgb 'violet' title ''
