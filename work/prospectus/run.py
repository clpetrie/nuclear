import os
import platform
import sys

os.system('pdflatex prospectus.tex')
if len(sys.argv) > 1:
   os.system('bibtex prospectus')
   os.system('pdflatex prospectus.tex')
   os.system('pdflatex prospectus.tex')
if platform.system() == 'Darwin':
   os.system('open prospectus.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open prospectus.pdf')
