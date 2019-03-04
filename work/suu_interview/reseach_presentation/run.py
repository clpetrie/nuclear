import os
import platform
import sys 

os.system('pdflatex research.tex')
if len(sys.argv) > 1:
   os.system('bibtex research')
   os.system('pdflatex research.tex')
   os.system('pdflatex research.tex')
if platform.system() == 'Darwin':
   os.system('open research.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open research.pdf')
os.system('rm *.bib *aux *key *log *nav *out *xml *snm *toc *bbl *blg')
