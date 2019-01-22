import os
import platform
import sys

os.system('pdflatex pd_interview.tex')
if len(sys.argv) > 1:
   os.system('bibtex pd_interview')
   os.system('pdflatex pd_interview.tex')
   os.system('pdflatex pd_interview.tex')
if platform.system() == 'Darwin':
   os.system('open pd_interview.pdf')
if platform.system() == 'Linux':
   os.system('gnome-open pd_interview.pdf')
