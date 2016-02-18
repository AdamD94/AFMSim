#locations plot

  
unset palette

set size ratio -1
set origin -0.033, 0

unset cblabel
set xlabel 'x [Å]'
set ylabel 'y [Å]'
set zlabel 'z [Å]'

splot 'Locations.dat' using 1:2:3:3 notitle palette pt 7