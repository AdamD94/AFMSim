#xy plot with potential 

set xlabel 'x [Å]'
set ylabel 'y [Å]'
set zlabel 'z [Å]'

set size ratio -1

set palette model RGB
set palette model RGB defined (	0 "dark-blue",2 "green",3 "yellow",4 "red",5 "black")

set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'F [nN]'

set size ratio -1

set parametric;
plot 'Surface.dat' using 1:2:4 notitle palette, 'Locations.dat' using 1:2:3 notitle pt 7