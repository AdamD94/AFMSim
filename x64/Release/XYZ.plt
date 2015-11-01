#xyz plot with layers coloured differently 

set palette model RGB defined (0 "red",1 "blue", 2 "green")
unset colorbox
set origin -0.4, -0.4                                     
set size 1.8, 1.8
set view equal xyz
set xlabel 'x [nm]'
set ylabel 'y [nm]'
set zlabel 'z [nm]'
set view 0,0
splot 'AFM.dat' using 1:2:3:3 notitle with points pt 7 palette