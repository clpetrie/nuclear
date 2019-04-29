import os
import platform
import sys

if(len(sys.argv) == 1):
   os.system('pdflatex dis.tex')
   os.system('bibtex dis')
   os.system('pdflatex dis.tex')
   os.system('pdflatex dis.tex')
else:
   os.system('pdflatex dis.tex')
if platform.system() == 'Darwin':
   os.system('open dis.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open dis.pdf')
