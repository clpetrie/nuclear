import os
import platform
import sys
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def main():
   #set up parameters
   num=40 #number of parameters
   filename='4he.vmc10000.opt1.out'

   if len(sys.argv) > 1:
      filename=str(sys.argv[1])
   #read in parameter data
   os.system('grep "new parameters" '+filename+' > temp.txt')
   params=np.loadtxt('temp.txt',usecols=(range(3,num+3)))
   os.system('rm temp.txt')
   #read in energy data
   os.system('grep "total energy" '+filename+' > temp.txt')
   energy=np.loadtxt('temp.txt',usecols=range(3,4))
   eerr=np.loadtxt('temp.txt',usecols=(range(5,6)))
   os.system('rm temp.txt')

   #do plot of energy then parameters
   axes = AxesSequence()
   for i, ax in zip(range(num+1), axes):
      if(i==0):
         x=range(1,len(energy)+1)
         ax.errorbar(x,energy,yerr=eerr,fmt='-o')
         ax.set_title('Energy')
      else:
         ax.plot(params[:,i-1])
         ax.set_title('Parameter {}'.format(i))
   axes.show()

class AxesSequence(object):
    """Creates a series of axes in a figure where only one is displayed at any
    given time. Which plot is displayed is controlled by the arrow keys."""
    def __init__(self):
        self.fig = plt.figure()
        self.axes = []
        self._i = 0 # Currently displayed axes index
        self._n = 0 # Last created axes index
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

    def __iter__(self):
        while True:
            yield self.new()

    def new(self):
        # The label needs to be specified so that a new axes will be created
        # instead of "add_axes" just returning the original one.
        ax = self.fig.add_axes([0.15, 0.1, 0.8, 0.8],
                               visible=False, label=self._n)
        self._n += 1
        self.axes.append(ax)
        return ax

    def on_keypress(self, event):
        if event.key == 'right':
            self.next_plot()
        elif event.key == 'left':
            self.prev_plot()
        else:
            return
        self.fig.canvas.draw()

    def next_plot(self):
        if self._i < len(self.axes):
            self.axes[self._i].set_visible(False)
            self.axes[self._i+1].set_visible(True)
            self._i += 1

    def prev_plot(self):
        if self._i > 0:
            self.axes[self._i].set_visible(False)
            self.axes[self._i-1].set_visible(True)
            self._i -= 1

    def show(self):
        self.axes[0].set_visible(True)
        plt.show()

if __name__ == '__main__':
    main()
