#set terminal pngcairo enhanced font 'Verdana,11'
set terminal png enhanced font 'Verdana,11'

sys=1 #1=alpha, 2=14n, 3=14n2p, 4=all_together
datafile='data.txt'
if (sys==1) imin="3"
if (sys==1) imax="8"
if (sys==2) imin="11"
if (sys==2) imax="16"
if (sys==3) imin="19"
if (sys==3) imax="24"

if (sys==1) set output 'alpha.png'
if (sys==2) set output '14n.png'
if (sys==3) set output '14n2p.png'
if (sys==4) set output 'plot.png'

lsize=2
psize=1.2

set xlabel "Density"
set ylabel "Alpha Particle Energy (MeV)"

if (sys<=3) \
   plot "<(sed -n '".imin.",".imax."p' data.txt)" u 1:2:3 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'Linear', \
      "<(sed -n '".imin.",".imax."p' data.txt)" u 1:2 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
      "<(sed -n '".imin.",".imax."p' data.txt)" u 1:6:7 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'IP Quadratic', \
      "<(sed -n '".imin.",".imax."p' data.txt)" u 1:6 with lines linetype 2 linewidth lsize lc rgb "blue" notitle

#      "<(sed -n '".imin.",".imax."p' data.txt)" u 1:4:5 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'IP Quadratic', \
#      "<(sed -n '".imin.",".imax."p' data.txt)" u 1:4 with lines linetype 2 linewidth lsize lc rgb "green" notitle, \
if (sys==4) \
   plot "<(sed -n '11,16p' data.txt)" u 1:2:3 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'Linear - 14n', \
      "<(sed -n '11,16p' data.txt)" u 1:2 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
      "<(sed -n '11,16p' data.txt)" u 1:4:5 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'IP - 14n', \
      "<(sed -n '11,16p' data.txt)" u 1:4 with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
      "<(sed -n '11,16p' data.txt)" u 1:6:7 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'IP (IP Opt) - 14n', \
      "<(sed -n '11,16p' data.txt)" u 1:6 with lines linetype 2 linewidth lsize lc rgb "green" notitle, \
      "<(sed -n '19,24p' data.txt)" u 1:2:3 with yerrorbars linetype 1 linewidth 1.2 pointtype 9 pointsize psize lc rgb "red" title 'Linear - 14n2p', \
      "<(sed -n '19,24p' data.txt)" u 1:2 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
      "<(sed -n '19,24p' data.txt)" u 1:4:5 with yerrorbars linetype 1 linewidth 1.2 pointtype 9 pointsize psize lc rgb "blue" title 'IP - 14n2p', \
      "<(sed -n '19,24p' data.txt)" u 1:4 with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
      "<(sed -n '19,24p' data.txt)" u 1:6:7 with yerrorbars linetype 1 linewidth 1.2 pointtype 9 pointsize psize lc rgb "green" title 'IP (IP Opt) - 14n2p', \
      "<(sed -n '19,24p' data.txt)" u 1:6 with lines linetype 2 linewidth lsize lc rgb "green" notitle

#   plot "<(sed -n '3,8p' data.txt)" u 1:2:3 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'Linear - alpha', \
#      "<(sed -n '3,8p' data.txt)" u 1:2 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
#      "<(sed -n '3,8p' data.txt)" u 1:4:5 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'IP - alpha', \
#      "<(sed -n '3,8p' data.txt)" u 1:4 with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
#      "<(sed -n '3,8p' data.txt)" u 1:6:7 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'IP (IP Opt) - alpha', \
#      "<(sed -n '3,8p' data.txt)" u 1:6 with lines linetype 2 linewidth lsize lc rgb "green" notitle, \
