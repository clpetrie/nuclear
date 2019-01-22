import os
import platform
import sys

os.system('pdflatex highlevel.tex')
if len(sys.argv) > 1:
   os.system('bibtex highlevel')
   os.system('pdflatex highlevel.tex')
   os.system('pdflatex highlevel.tex')
if platform.system() == 'Darwin':
   os.system('open highlevel.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open highlevel.pdf')
