import os

#define file names
infile='data.out'
outfile='data.dat'

#make outfile
os.system("grep 'energy for block =' "+infile+" > "+outfile)
os.system("gnuplot plot.gp")
os.system("open data.gif")
