import os
import platform

os.system('pdflatex quadcorr.tex')
os.system('bibtex quadcorr')
os.system('pdflatex quadcorr.tex')
os.system('pdflatex quadcorr.tex')
if platform.system() == 'Darwin':
   os.system('open quadcorr.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open quadcorr.pdf')
