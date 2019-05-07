import os
import platform
import sys

if(len(sys.argv) == 1):
   os.system('pdflatex methods.tex')
   os.system('bibtex methods')
   os.system('pdflatex methods.tex')
   os.system('pdflatex methods.tex')
else:
   os.system('pdflatex methods.tex')
if platform.system() == 'Darwin':
   os.system('open methods.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open methods.pdf')
