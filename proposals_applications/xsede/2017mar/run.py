import os
import platform

os.system('pdflatex xsede.tex')
os.system('bibtex xsede')
os.system('pdflatex xsede.tex')
os.system('pdflatex xsede.tex')
if platform.system() == 'Darwin':
   os.system('open xsede.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open xsede.pdf')
