set title "Probability distribution for n=4 and n =6"
set xlabel "X" #font ", 16"
set ylabel "Occurance in The interval" #font ", 16"
set grid
set xr[-2:2]
plot '4_sum' w lp, '6_sum' w lp
set terminal png
set output "blas3.png"
replot
