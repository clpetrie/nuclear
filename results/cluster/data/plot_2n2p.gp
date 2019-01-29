#set terminal pngcairo enhanced font 'Verdana,11'
set terminal png enhanced font 'Verdana,11'

datafile='data.txt'
#datafile='vke.txt'
#datafile='vc.txt'
#datafile='vt.txt'
#datafile='vs.txt'
#datafile='vst.txt'
#datafile='vten.txt'
#datafile='vtent.txt'
imin="4"
imax="14"

set output 'alpha_and_2n2p.png'

lsize=2
psize=1.2

set xlabel "Density"
set ylabel "Energy (MeV)"

plot "<(sed -n '".imin.",".imax."p' ".datafile.")" u 1:2:3 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'Linear - Cluster', \
   "<(sed -n '".imin.",".imax."p' ".datafile.")" u 1:2 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
   "<(sed -n '".imin.",".imax."p' ".datafile.")" u 1:6:7 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'IP - Cluster', \
   "<(sed -n '".imin.",".imax."p' ".datafile.")" u 1:6 with lines linetype 2 linewidth lsize lc rgb "green" notitle, \
   "<(sed -n '".imin.",".imax."p' ../alpha/data/".datafile.")" u 1:2:3 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "orange" title 'Linear - 2n2p', \
   "<(sed -n '".imin.",".imax."p' ../alpha/data/".datafile.")" u 1:2 with lines linetype 2 linewidth lsize lc rgb "orange" notitle, \
   "<(sed -n '".imin.",".imax."p' ../alpha/data/".datafile.")" u 1:4:5 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'IP - 2n2p', \
   "<(sed -n '".imin.",".imax."p' ../alpha/data/".datafile.")" u 1:4 with lines linetype 2 linewidth lsize lc rgb "blue" notitle
