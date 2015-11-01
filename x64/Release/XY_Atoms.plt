#xy plot with basis atoms coloured differently to lattice atoms

set palette model RGB defined (0 "red",1 "blue", 2 "green")
unset colorbox
unset size
unset origin
set size ratio -1  
set xlabel 'x [nm]'
set ylabel 'y [nm]'
plot 'AFM.dat' using 1:2:6 notitle with points pt 7 palette
