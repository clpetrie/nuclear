import os
import platform

os.system('pdflatex paper.tex')
os.system('bibtex paper')
os.system('pdflatex paper.tex')
os.system('pdflatex paper.tex')
if platform.system() == 'Darwin':
   os.system('open paper.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open paper.pdf')
