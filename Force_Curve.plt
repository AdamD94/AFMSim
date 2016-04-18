#xy plot with potential 
set term pngcairo enhanced color size 1024,798 crop font 'Veranda, 18'
set output 'Force_Curve.png'
  
set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "blue",4 "light-blue",5 "light-green",6 "yellow", 8 "red", 10"black")

set size ratio -5
set origin -0.033, 0

set xtics
set ytics

unset cblabel
set xlabel 'x [Å]'
set ylabel 'z [Å]'
set cblabel 'F [nN]'

set pm3d map interpolate 2,2
splot 'Force_Curve.dat' using 1:3:4 notitle
