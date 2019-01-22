import os
import platform
import sys

os.system('pdflatex int.tex')
if len(sys.argv) > 1:
   os.system('bibtex int')
   os.system('pdflatex int.tex')
   os.system('pdflatex int.tex')
if platform.system() == 'Darwin':
   os.system('open int.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open int.pdf')
