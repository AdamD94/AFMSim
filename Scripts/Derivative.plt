#Tip Height Plot
unset xrange
unset yrange
unset zrange
unset cblabel


set term pngcairo enhanced color size 1024,798 crop font 'Veranda, 18'
set output 'Derivative.png'
  
set palette model RGB
set palette rgb 34,35,36; 
set palette functions gray, gray, gray

set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "blue",4 "light-blue",5 "light-green",6 "yellow", 8 "red", 10"black")

set size ratio -1
set origin -0.033, 0

unset cbtics
unset colorbox
set cblabel 'z [Å]'

set lmargin 0
set tmargin 0
set bmargin 0

unset xtics
unset ytics
unset border


set pm3d map interpolate 2,2
splot 'Topology.dat' using 1:2:4:4 notitle palette