set terminal postscript eps color 20 lw 3
set output '| epstopdf -f -o=computation-time.pdf'

set key autotitle columnheader

set ylabel "Time (sec)"
set xlabel "Path lengths"
set key left

# set xtics (0, '' 1, 2, 4, 8, 16, 32, 64)
# set ytics (4, 16, 64, 256, 1024, 4096, "16,384" 16384)
# set logscale y 2

set xtics ('5,000,000' 5000000, '10,000,000' 10000000, '15,000,000' 15000000, '20,000,000' 20000000, '25,000,000' 25000000)

# set auto x

f(x) = m*x + b
fit f(x) 'computation-time-path-homotopy.dat' using 3:5 via m,b

plot 'computation-time-path-homotopy.dat' using 3:5 with points title "Homotopy test", f(x) title 'Model Fit'

