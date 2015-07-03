#!/bin/bash
# int-or-string.sh

arg2=0
arg3=1

for i in `seq 1 20`; do

arg1=$i

msub -N toffStep${arg1}coef${arg2}${arg3} -q normal-long python OptimizationToffoli.py ${arg1} ${arg2} ${arg3} 

done
