#!/bin/bash
# int-or-string.sh

for i in `seq 10 99`; do 

echo "set term pdf
set output 'out${i}.pdf'
plot 'toffoli_optimized0.00${i}coef0.50.5' w l" > plot
gnuplot plot 

done
