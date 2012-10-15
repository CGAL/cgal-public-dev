set title "Effect of Randomization"
set xlabel "No. of points"
set ylabel "Time (sec)"

plot "randomization_influence_on_points.plot" u 1:2 t 'Randomized' w l lw 5, \
"randomization_influence_on_sorted_points.plot" u 1:2 t 'Deteministic' w l lw 5

load 'save.plt'