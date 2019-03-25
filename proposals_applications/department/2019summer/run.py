import os
import platform

os.system('pdflatex petrie_proposal.tex')
if platform.system() == 'Darwin':
   os.system('open petrie_proposal.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open petrie_proposal.pdf')
