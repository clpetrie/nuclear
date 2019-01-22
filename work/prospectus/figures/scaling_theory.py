import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv

fs=20
lw=5
xmin=1
xmax=28
A=np.arange(xmin,xmax)
lin=A*(A-1)/2.0
ip=A*(A-1)*(A-2)*(A-3)/8.0
quad=lin**2-lin
quadfix=quad-ip

yip=(lin+ip)/lin
yquad=(lin+quad)/lin
yquadfix=(lin+quadfix)/lin

plt.ylabel('Number Correlation Terms',fontsize=fs)
plt.xlabel('A',fontsize=fs)
plt.xlim([xmin,xmax])

#plt.plot(A,ylin,'b',marker='-',A,yip,'r.-',A,yquad,'g:',linewidth=lw)
plt.plot(A,yip,'b',A,yquad,'r',A,yquadfix,'g',linewidth=lw)
plt.figure()
plt.plot(A,yquad/yip,'r',A,yquadfix/yip,'g')

plt.savefig('scaleing_theory.pdf',format='pdf')
plt.show()
