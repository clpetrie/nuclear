import os
import platform

os.system('pdflatex research.tex')
if platform.system() == 'Darwin':
   os.system('open research.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open research.pdf')
os.system('rm *.bib *aux *key *log *nav *out *xml *snm *toc')
