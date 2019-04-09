import os
import platform
import sys

os.system('pdflatex april_petrie.tex')
if len(sys.argv) > 1:
   os.system('bibtex april_petrie')
   os.system('pdflatex april_petrie.tex')
   os.system('pdflatex april_petrie.tex')
if platform.system() == 'Darwin':
   os.system('open april_petrie.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open april_petrie.pdf')
