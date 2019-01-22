#!/bin/bash
min=1
max=25
for ((i=$min;i<=$max;i++));
do
   grep 'total energy' he4lin.unc10000_$i.out > dat/he4lin.unc10000_$i.dat
done
