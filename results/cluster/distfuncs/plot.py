import os
import platform
import sys
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

file0005ip='temp0005ip.dmc'
file001ip='temp001ip.dmc'
file002ip='temp002ip.dmc'
file003ip='temp003ip.dmc'

r0005ip=np.loadtxt(file0005ip,usecols=range(0,1))
pp0005ip=np.loadtxt(file0005ip,usecols=range(5,6))
dpp0005ip=np.loadtxt(file0005ip,usecols=range(6,7))
nn0005ip=np.loadtxt(file0005ip,usecols=range(7,8))
dnn0005ip=np.loadtxt(file0005ip,usecols=range(8,9))

r001ip=np.loadtxt(file001ip,usecols=range(0,1))
pp001ip=np.loadtxt(file001ip,usecols=range(5,6))
dpp001ip=np.loadtxt(file001ip,usecols=range(6,7))
nn001ip=np.loadtxt(file001ip,usecols=range(7,8))
dnn001ip=np.loadtxt(file001ip,usecols=range(8,9))

r002ip=np.loadtxt(file002ip,usecols=range(0,1))
pp002ip=np.loadtxt(file002ip,usecols=range(5,6))
dpp002ip=np.loadtxt(file002ip,usecols=range(6,7))
nn002ip=np.loadtxt(file002ip,usecols=range(7,8))
dnn002ip=np.loadtxt(file002ip,usecols=range(8,9))

r003ip=np.loadtxt(file003ip,usecols=range(0,1))
pp003ip=np.loadtxt(file003ip,usecols=range(5,6))
dpp003ip=np.loadtxt(file003ip,usecols=range(6,7))
nn003ip=np.loadtxt(file003ip,usecols=range(7,8))
dnn003ip=np.loadtxt(file003ip,usecols=range(8,9))

#plt.plot(r0005ip[0:50],pp0005ip[0:50],'ro',r003ip[0:50],pp003ip[0:50],'bo')
plt.xlim((5,35))
plt.ylim((0,0.05))
plt.plot(r0005ip,pp0005ip,'ro',r001ip,pp001ip,'go',r002ip,pp002ip,'yo',r003ip,pp003ip,'bo')
#plt.plot(r0005ip,pp0005ip,'r-',r001ip,pp001ip,'g-',r002ip,pp002ip,'y-',r003ip,pp003ip,'b-')
#plt.plot(r0005ip,nn0005ip,'ro',r001ip,nn001ip,'go',r002ip,nn002ip,'yo',r003ip,nn003ip,'bo')
plt.show()
