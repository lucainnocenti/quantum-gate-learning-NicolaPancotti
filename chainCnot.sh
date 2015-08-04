#!/bin/bash
# int-or-string.sh

arg2=0.5
arg3=0.5

for i in `seq 5 100`; do

arg1=$i

msub -N cnot${arg1}coef${arg2}${arg3} -q normal python OptimizationChainHeis.py ${arg1} ${arg2} ${arg3} 7

done
