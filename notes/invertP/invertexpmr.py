import matplotlib.pyplot as plt
import numpy as np
import random as ran

num=10000;
ranr=[]
r=[]

for i in range(0,num):
   temp=np.log(1/ran.random())
   ranr.append(temp)

P,edges=np.histogram(ranr,bins=50)
P=P/(float(np.sum(P))*(edges[1]-edges[0])) #normalize P to add to 1
for i in range(0,len(edges)-1):
   r.append(0.5*(edges[i+1]+edges[i]))
plt.plot(r,P)
r=np.arange(0,10,0.1)
plt.plot(r,np.exp(-r))
plt.show()
