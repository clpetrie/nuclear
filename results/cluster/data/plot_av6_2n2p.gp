#set terminal pngcairo enhanced font 'Verdana,11'
set terminal png enhanced font 'Verdana,11'

corr=2 #1=linear, 2=ip
imin="4"
imax="14"
set output 'av6_2n2p.png'
if (corr==1) \
   cmin1=2; \
   cmax1=3
if (corr==2) \
   cmin1=4; \
   cmax1=5

fileenergy='../alpha/data/data.txt'
fileke='../alpha/data/vke.txt'
filec='../alpha/data/vc.txt'
filet='../alpha/data/vt.txt'
files='../alpha/data/vs.txt'
filest='../alpha/data/vst.txt'
fileten='../alpha/data/vten.txt'
filetent='../alpha/data/vtent.txt'

lsize=2
psize=1.2

set xlabel "Density"
set ylabel "Energy (MeV)"
set xrange[0:0.0135]

plot "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'Etotal', \
   "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'vc', \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "orange" title 'vt', \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "orange" notitle, \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "violet" title 'vs', \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "violet" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "cyan" title 'vst', \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "cyan" notitle, \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "black" title 'vten', \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "black" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'vtent', \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "green" notitle

#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "yellow" title 'vke', \
#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "yellow" notitle, \
