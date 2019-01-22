import os
import platform
import sys

os.system('pdflatex talk.tex')
if len(sys.argv) > 1:
   os.system('bibtex talk')
   os.system('pdflatex talk.tex')
   os.system('pdflatex talk.tex')
if platform.system() == 'Darwin':
   os.system('open talk.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open talk.pdf')
