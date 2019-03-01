import os
import platform

os.system('pdflatex radiation.tex')
if platform.system() == 'Darwin':
   os.system('open radiation.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open radiation.pdf')
os.system('rm *.bib *aux *key *log *nav *out *xml *snm *toc')
