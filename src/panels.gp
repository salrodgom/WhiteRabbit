set terminal postscript eps color enhanced 'Helvetica,12'
set output 'panels.eps'
set encoding iso_8859_1
# size
set size 0.5,1
set multiplot
set xrange [0.0:*]
# First plot  (large)
set origin 0,0
set size 0.5,0.5
unset key
set xlabel "{/Symbol d}^{(8-ring)} [{/E \305}]"
unset ytics
#set xtics 1
set xtics -0.0,0.2,2.0
plot "8-histogram.txt" u 1:2 w l lt 1 lw 1.5 lc rgb"#FF6347"
####################################################################### 6,4
set origin 0.0,0.5
set size 0.25,0.5
set xtics -0.0,0.2,1.6
set xlabel "{/Symbol d}^{(4-ring)} [{/E \305}]"
plot "4-histogram.txt" u 1:2 w l lt 1 lw 1.5 lc rgb"#FF6347"
######
set origin 0.25,0.5
set size 0.25,0.5
set xtics -0.0,0.2,1.6
set xlabel "{/Symbol d}^{(6-ring)} [{/E \305}]"
plot "6-histogram.txt" u 1:2:(1500) w l lt 1 lw 1.5 lc rgb"#FF6347"
unset multiplot
