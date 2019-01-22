import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv

fs=20
lw=3
xmin=-1
xmax=1
x=np.arange(xmin,xmax,0.05)
#y=1/(np.sqrt(2*3.14159))*np.exp(-x**2/2)
#y=-np.sqrt(2)*erfinv(2*x)
y=erfinv(2*x)

plt.ylabel('P(r)',fontsize=fs)
plt.xlabel('r',fontsize=fs)
plt.text(1.0, 0.5, r'$P(r)=e^{-r}$', fontsize=fs)
plt.xlim([xmin,xmax])

plt.plot(x,y,'b',linewidth=lw)

rand=np.random.rand(100000)
print -np.sqrt(2)*erfinv(2*0.5)
n, bins, patches = plt.hist(-np.sqrt(2)*erfinv(2*rand), 50, normed=1, facecolor='r')
plt.savefig('gausscdf.pdf',format='pdf')
plt.show()
