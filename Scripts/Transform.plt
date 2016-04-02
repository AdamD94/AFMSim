#Tip Height Plot
unset xrange
unset yrange

set term pngcairo enhanced color size 1024,798 crop font 'Veranda, 18'
set output 'Transform.png'
set xrange[0:4E8]

plot 'Transform.dat' with lines notitle