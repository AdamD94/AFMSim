#locations plot

unset palette

set size ratio -1

unset cblabel
set xlabel 'x [Å]'
set ylabel 'y [Å]'
set zlabel 'z [Å]'

plot 'Locations.dat' using 1:2 notitle pt 7