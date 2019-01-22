import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt

if sys.argv[1] != "":
   filename=[]
   num=len(sys.argv)-1
   for n in range(1,num+1):
      filename.append(sys.argv[n])

for n in range(0,num):
   data=np.loadtxt(filename[n],usecols=(3,))
   error=np.loadtxt(filename[n],usecols=(5,))

   print 'E = ',np.mean(data),' +/- ',np.sqrt(np.mean(np.square(error))),' MeV, ',filename[n]

plt.plot(data)
plt.show()
