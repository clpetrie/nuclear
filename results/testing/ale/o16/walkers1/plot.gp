set term x11 persist

doval='all'
file1='o16.ale_sites'
file2='o16.t2ale.alpha100_sites'
file3='o16.t2ale.alpha10000_sites'
file4='o16.t2ale.alpha1E32_sites'
post='.t2'
mylw=2

set logscale y
set format y "%E"
#plot file.'.t2' using (abs($6)) with lines, file.'.tz' using (abs($6)) with lines, file.'.j2' using (abs($6)) with lines, file.'.jz' using (abs($6)) with lines


plot file1.post using (abs($6)) with lines title 'no T2' lw mylw, file2.post using (abs($6)) with lines title 'alpha=100' lw mylw, file3.post using (abs($6)) with lines title 'alpha=10000' lw mylw, file4.post using (abs($6)) with lines title 'alpha=1E32' lw mylw
