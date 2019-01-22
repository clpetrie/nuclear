import os
import platform

#run latex and bibtex on my thesis
os.system('pdflatex expcorr.tex')
#os.system('bibtex expcorr')
#os.system('pdflatex expcorr.tex')
#os.system('pdflatex expcorr.tex')
if platform.system() == 'Darwin':
   os.system('open expcorr.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open expcorr.pdf')
