import os
import platform

os.system('pdflatex petrie.tex')
os.system('bibtex petrie')
os.system('pdflatex petrie.tex')
os.system('pdflatex petrie.tex')
if platform.system() == 'Darwin':
   os.system('open petrie.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open petrie.pdf')
