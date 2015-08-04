#!/bin/bash
# int-or-string.sh

arg2=0.5
arg3=0.5

for i in `seq 1 100`; do

arg1=$i

msub -N toffStep${arg1}coef${arg2}${arg3} -q normal python OptimizationFred.py ${arg1} ${arg2} ${arg3} 

done
