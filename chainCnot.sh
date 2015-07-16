#!/bin/bash
# int-or-string.sh

arg2=0.5
arg3=0.5

for i in `seq 5 50`; do

arg1=$i

msub -N toff${arg1}coef${arg2}${arg3} -q normal-long python OptimizationToffoli.py ${arg1} ${arg2} ${arg3} 

done
