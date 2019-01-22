import os
import platform

os.system('pdflatex outline.tex')
os.system('bibtex outline')
os.system('pdflatex outline.tex')
os.system('pdflatex outline.tex')
if platform.system() == 'Darwin':
   os.system('open outline.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open outline.pdf')
