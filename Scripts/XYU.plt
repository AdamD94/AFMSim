#xy plot with potential 
unset xrange
unset yrange
unset zrange
unset cbrange


set term pngcairo enhanced color size 1024,796 crop font 'Veranda, 18'
set output 'Surface.png'
  
set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "blue",4 "light-blue",5 "light-green",6 "yellow", 8 "red", 10"black")

set size ratio -1


unset cblabel
set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'F [nN]'


set lmargin 0
set tmargin 0
set bmargin 0

unset xtics
unset ytics
unset border


set pm3d map interpolate 2,2
splot 'Surface.dat' using 1:2:4:4 notitle palette