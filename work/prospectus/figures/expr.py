import numpy as np
import matplotlib.pyplot as plt

fs=20
lw=3
x=np.arange(0,10,0.05)
y=np.exp(-x)

plt.ylabel('P(r)',fontsize=fs)
plt.xlabel('r',fontsize=fs)
plt.text(1.0, 0.5, r'$P(r)=e^{-r}$', fontsize=fs)
plt.xlim([0,10])

plt.plot(x,y,'b',linewidth=lw)

rand=np.random.rand(100000)
n, bins, patches = plt.hist(np.log(1/rand), 50, normed=1, facecolor='r')
plt.savefig('expr.pdf',format='pdf')
plt.show()
