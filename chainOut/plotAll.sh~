#!/bin/bash
# int-or-string.sh

for i in `seq 10 50`; do 

echo "set term pdf
set output 'out${i}.pdf'
plot 'cnotHeis0.0${i}coef0.50.5N7' w l" > plot
gnuplot plot 

done

