import os
import platform

#run latex and bibtex on my thesis
os.system('pdflatex notes.tex')
os.system('bibtex notes')
os.system('pdflatex notes.tex')
os.system('pdflatex notes.tex')
if platform.system() == 'Darwin':
   os.system('open notes.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open notes.pdf')
