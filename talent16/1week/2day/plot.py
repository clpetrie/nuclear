import os

#define file names
infile='data.out'
outfile='energy.dat'

#make outfile
os.system("grep 'energy for block =' "+infile+" > "+outfile)
os.system("gnuplot plot.gp")
os.system("open energy.gif")
