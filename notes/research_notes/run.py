import os
import platform
import sys

os.system('pdflatex notes.tex')
if len(sys.argv) > 1:
   os.system('bibtex notes')
   os.system('pdflatex notes.tex')
   os.system('pdflatex notes.tex')
if platform.system() == 'Darwin':
   os.system('open notes.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open notes.pdf')
