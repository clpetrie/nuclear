from functions import *
import os
#import matplotlib.pyplot as plt
from matplotlib import pylab
import sys

#Define variables
searchfor="'total energy'"
imin=1
if len(sys.argv) > 1:
   filename=sys.argv[1]
else:
   filename='../he4/he4lin.dmc10000.out'

os.system("grep "+searchfor+" "+filename+" > temp.dat")
#Get data and determine uncertainties for blocks
mdata=[]
edata=[]
f = open("temp.dat")
for line in f.readlines():
   mdata.append(float(line[47:64]))
   edata.append(float(line[71:-4]))
f.close()
os.system("rm temp.dat")

pylab.xlabel('block')
pylab.ylabel('energy')
pylab.plot(mdata,'bo-')
pylab.show()
