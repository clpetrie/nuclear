#set terminal pngcairo enhanced font 'Verdana,11'
#set terminal png enhanced font 'Verdana,11' size 640,480
set terminal png enhanced font 'Verdana,11' size 740,480

corr1=1 #1=linear, 2=ip (lin opt), 3=ip (ip opt)
corr2=3 #1=linear, 2=ip (lin opt), 3=ip (ip opt)
sys=1 #1=alpha, 2=14n, 3=14n2p
if (sys==1) \
   imin="4"; \
   imax="15"; \
   set output 'av6_alpha.png'
if (sys==2) \
   imin="18"; \
   imax="29"; \
   set output 'av6_14n.png'
if (sys==3) \
   imin="32"; \
   imax="43"; \
   set output 'av6_14n2p.png'
if (corr1==1) \
   cmin1=2; \
   cmax1=3
if (corr1==2) \
   cmin1=4; \
   cmax1=5
if (corr1==3) \
   cmin1=6; \
   cmax1=7
if (corr2==1) \
   cmin2=2; \
   cmax2=3
if (corr2==2) \
   cmin2=4; \
   cmax2=5
if (corr2==3) \
   cmin2=6; \
   cmax2=7


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

set xlabel '{/Symbol r} (fm^{-3})' enhanced
set ylabel 'E (MeV)' enhanced
#set xrange[0:0.0140]

#set xrange[0:0.01]
#set size 0.78,0.78
#set key at 0.014,18

set xrange[0:0.01]
set size 0.85,1.00
set key at 0.0128,18

plot "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "red" title 'E_{total} - lin', \
   "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "blue" title 'v_{c} - lin', \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "orange" title 'v_{{/Symbol t}} - lin', \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "orange" notitle, \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "violet" title 'v_{{/Symbol s}} - lin', \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "violet" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "cyan" title 'v_{{/Symbol s}{/Symbol t}} - lin', \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "cyan" notitle, \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "black" title 'v_{S} - lin', \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "black" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "green" title 'v_{S{/Symbol t}} - lin', \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "green" notitle, \
   "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "red" title 'E_{total} - ip', \
   "<(sed -n '".imin.",".imax."p' ".fileenergy.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "red" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "blue" title 'v_{c} - ip', \
   "<(sed -n '".imin.",".imax."p' ".filec.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "blue" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "orange" title 'v_{{/Symbol t}} - ip', \
   "<(sed -n '".imin.",".imax."p' ".filet.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "orange" notitle, \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "violet" title 'v_{{/Symbol s}} - ip', \
   "<(sed -n '".imin.",".imax."p' ".files.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "violet" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "cyan" title 'v_{{/Symbol s}{/Symbol t}} - ip', \
   "<(sed -n '".imin.",".imax."p' ".filest.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "cyan" notitle, \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "black" title 'v_{S} - ip', \
   "<(sed -n '".imin.",".imax."p' ".fileten.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "black" notitle, \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 3 pointsize psize lc rgb "green" title 'v_{S{/Symbol t}} - ip', \
   "<(sed -n '".imin.",".imax."p' ".filetent.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "green" notitle

#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin1:cmax1 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "yellow" title 'vke', \
#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin1 with lines linetype 2 linewidth lsize lc rgb "yellow" notitle, \
#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin2:cmax2 with yerrorbars linetype 1 linewidth 1.2 pointtype 7 pointsize psize lc rgb "yellow" title 'vke', \
#   "<(sed -n '".imin.",".imax."p' ".fileke.")" u 1:cmin2 with lines linetype 2 linewidth lsize lc rgb "yellow" notitle, \
