import os
import platform

os.system('pdflatex diss.tex')
os.system('bibtex diss')
os.system('pdflatex diss.tex')
os.system('pdflatex diss.tex')
if platform.system() == 'Darwin':
   os.system('open diss.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open diss.pdf')
