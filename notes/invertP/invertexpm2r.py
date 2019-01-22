import matplotlib.pyplot as plt
import numpy as np
import random as ran
from math import pi

num=100000
rmin=0.0
rmax=10.0
a=0.0
b=0.5 #max value of P
nbin=50
ranr=[]
r=[]
mybin=[]

for i in range(0,num):
   temp=0.5*np.log(1/(2*ran.uniform(a,b))) #P=2*exp(-2r)
   ranr.append(temp)

step=(rmax-rmin)/nbin
for i in range(0,nbin+1):
   mybin.append(i*step)

P,edges=np.histogram(ranr,bins=mybin)
#P=P/(float(np.sum(P))*(edges[1]-edges[0])) #normalize P to add to 1
#for i in range(0,len(edges)-1):
#   r.append(0.5*(edges[i+1]+edges[i]))
P=P/(float(np.sum(P))*(mybin[1]-mybin[0])) #normalize P to add to 1
for i in range(0,len(mybin)-1):
   r.append(0.5*(mybin[i+1]+mybin[i]))
plt.plot(r,P,linewidth=3)
print r
r=np.arange(rmin,rmax,0.05)
plt.plot(r,2.0*np.exp(-2*r),'r',linewidth=2)
#plt.xlim(rmin,rmax)
plt.show()
