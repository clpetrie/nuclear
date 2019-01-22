import os
import platform
import sys

os.system('pdflatex 4c.tex')
if len(sys.argv) > 1:
   os.system('bibtex 4c')
   os.system('pdflatex 4c.tex')
   os.system('pdflatex 4c.tex')
if platform.system() == 'Darwin':
   os.system('open 4c.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open 4c.pdf')
