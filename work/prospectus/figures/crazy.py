import numpy as np
import matplotlib.pyplot as plt

fs=20
lw=3
x=np.arange(0,10,0.05)
y=(1.9*np.sin(2*x)+3*np.cos(0.5*x)-2.3*np.sin(3*x)-2.2*np.cos(x))**2/121.324 # crazy, norm 0 to 10

plt.ylabel('P(r)',fontsize=fs)
plt.xlabel('r',fontsize=fs)
plt.text(2.0, 0.25, r'$P(r)=\left|\Psi(r)\right|^2$', fontsize=fs)

plt.plot(x,y,'b',linewidth=lw)
plt.savefig('crazy.pdf',format='pdf')
plt.show()
