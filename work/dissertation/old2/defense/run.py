import os
import platform
import sys

os.system('pdflatex defense.tex')
if len(sys.argv) > 1:
   os.system('bibtex defense')
   os.system('pdflatex defense.tex')
   os.system('pdflatex defense.tex')
if platform.system() == 'Darwin':
   os.system('open defense.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open defense.pdf')
