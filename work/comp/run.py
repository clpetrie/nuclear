import os
import platform

os.system('pdflatex comp.tex')
os.system('bibtex comp')
os.system('pdflatex comp.tex')
os.system('pdflatex comp.tex')
if platform.system() == 'Darwin':
   os.system('open comp.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open comp.pdf')
