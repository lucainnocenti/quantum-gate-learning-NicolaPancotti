#!/bin/bash
# int-or-string.sh

for i in `seq 10 99`; do 

echo "set term pdf
set xlabel 'optimization steps'
set ylabel 'Likelihood Function'
set title 'stochastic gradient descent'
set yrange [0:1.1]
set output 'out${i}.pdf'
plot 'fred_optimized0.00${i}coef0.50.5' w l title 'optimization Fredkin gate'" > plot
gnuplot plot

done

rm plot
