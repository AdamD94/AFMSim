#xyz plot with potential, xy and xzplanes
set term pngcairo enhanced color size 1024,796 crop font 'Veranda, 18'
set output 'Combination.png'
  
set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "blue",4 "light-blue",5 "light-green",6 "yellow", 8 "red", 10"black")

set size ratio -1

set cblabel 'F [nN]'

unset xlabel
unset ylabel
unset zlabel
unset border
unset xtics
unset ytics
unset ztics
set hidden3d
splot 'Surface.dat' using 1:2:3:4 notitle palette, 'Force_Curve.dat' using 1:2:3:4 notitle palette
