#set terminal pngcairo enhanced font 'Verdana,11'
set terminal png enhanced font 'Verdana,11'

corr=3 #1=linear, 2=ip (lin opt), 3=ip (ip opt)
sys=3 #1=alpha, 2=14n, 3=14n2p
if (sys==1) \
   imin="4"; \
   imax="14"; \
   set output 'av6_alpha.png'
if (sys==2) \
   imin="17"; \
   imax="27"; \
   set output 'av6_14n.png'
if (sys==3) \
   imin="30"; \
   imax="40"; \
   set output 'av6_14n2p.png'
if (corr==1) \
   cmin=2; \
   cmax=3
if (corr==2) \
   cmin=4; \
   cmax=5
if (corr==3) \
   cmin=6; \
   cmax=7

fileenergy='data.txt'
fileke='vke.txt'
filec='vc.txt'
filet='vt.txt'
files='vs.txt'
filest='vst.txt'
fileten='vten.txt'
filetent='vtent.txt'

lsize=2
psize=1.2

set xlabel "Density"
set ylabel "Energy (MeV)"
set xrange[0:0.013]

plot "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'Etotal', \
   "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'vc', \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "orange" title 'vt', \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "orange" notitle, \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "violet" title 'vs', \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "violet" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "cyan" title 'vst', \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "cyan" notitle, \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "black" title 'vten', \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "black" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'vtent', \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "green" notitle
#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin:cmax with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "yellow" title 'vke', \
#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin with lines linetype 2 linewidth lsize lc rgb "yellow" notitle, \
