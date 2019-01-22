from functions import *
import os
import matplotlib.pyplot as plt
import sys

#Define variables
searchfor="'total energy'"
imin=1
if len(sys.argv) > 1:
   filename=sys.argv[1]
else:
   filename='unconstrained/he4lin_dmc10000.out'
N=100 #number of bunching intervals
fs=17
uncertainty=[]
version=1 #0=old, 1=new

os.system("grep "+searchfor+" "+filename+" > temp.dat")

#Get data and determine uncertainties for blocks
rawmean=[]
rawerr=[]
f = open("temp.dat")
for line in f.readlines():
   rawmean.append(float(line[47:64]))
   rawerr.append(float(line[71:-4]))
f.close()
mdata=rawmean[1:len(rawmean)]
edata=rawerr[1:len(rawerr)]
os.system("rm temp.dat")
