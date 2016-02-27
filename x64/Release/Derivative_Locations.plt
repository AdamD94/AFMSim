#xy plot with potential 

  
set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "blue",4 "light-blue",5 "light-green",6 "yellow", 8 "red", 10"black")

set size ratio -1
set origin -0.033, 0

unset cblabel
set xlabel 'x [Å]'
set ylabel 'z [Å]'
set cblabel 'dz'


set parametric;
plot 'Topology.dat' using 1:2:4:4 palette,  'VestaObj.dat' using 1:2:3 notitle pt 7
