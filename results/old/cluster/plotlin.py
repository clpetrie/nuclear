import os 
import platform


num=2
filename='lin.dmc'
prefix1='14n.'
prefix2='14n2p.'

#os.system('echo "set term pdf size 4,6" >> temp.gp')
os.system('echo "set term postscript size 4,6" >> temp.gp')
os.system('echo "set output \'lin.eps\'" >> temp.gp')
os.system('echo "set multiplot layout 3,2" >> temp.gp')
#os.system('echo "set yrange [-0.2:0.5]" >> temp.gp')
os.system('echo "set title \'rhopup\'" >> temp.gp')
os.system('echo "plot \''+prefix1+filename+'\' using 1:2 title \'14n\' with lines, \''+prefix2+filename+'\' using 1:2 with lines title \'14n2p\'" >> temp.gp')
os.system('echo "set title \'rhopdown\'" >> temp.gp')
os.system('echo "plot \''+prefix1+filename+'\' using 1:4 title \'14n\' with lines, \''+prefix2+filename+'\' using 1:4 with lines title \'14n2p\'" >> temp.gp')
os.system('echo "set title \'rhonup\'" >> temp.gp')
os.system('echo "plot \''+prefix1+filename+'\' using 1:6 title \'14n\' with lines, \''+prefix2+filename+'\' using 1:6 with lines title \'14n2p\'" >> temp.gp')
os.system('echo "set title \'rhondown\'" >> temp.gp')
os.system('echo "plot \''+prefix1+filename+'\' using 1:8 title \'14n\' with lines, \''+prefix2+filename+'\' using 1:8 with lines title \'14n2p\'" >> temp.gp')
os.system('echo "set title \'rhop\'" >> temp.gp')
os.system('echo "plot \''+prefix1+filename+'\' using 1:10 title \'14n\' with lines, \''+prefix2+filename+'\' using 1:10 with lines title \'14n2p\'" >> temp.gp')
os.system('echo "set title \'rhon\'" >> temp.gp')
os.system('echo "plot \''+prefix1+filename+'\' using 1:12 title \'14n\' with lines, \''+prefix2+filename+'\' using 1:12 with lines title \'14n2p\'" >> temp.gp')
os.system('echo "unset multiplot" >> temp.gp')
os.system('gnuplot temp.gp')

os.system('rm temp.gp')
if platform.system() == 'Darwin':
   os.system('open lin.eps')
if platform.system() == 'Linux':
   os.system('gnome-open lin.eps')
