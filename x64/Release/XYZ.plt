#xy plot with potential 

set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "blue",4 "light-blue",5 "light-green",6 "yellow", 8 "red", 10"black")

set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'F [nN]'

set size ratio -1

set parametric;
plot 'Surface.dat' using 1:2:4 notitle palette, 'Locations.dat' using 1:2:3 notitle pt 7 font "20" 