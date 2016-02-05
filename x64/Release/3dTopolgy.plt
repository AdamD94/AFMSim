#Tip Height Plot

set palette model RGB
set palette rgb 34,35,36; 

set view equal xyz
set origin -0.033, 0

unset cblabel
set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'z [Å]'

set pm3d interpolate 2,2
splot 'Tip_Height.dat' using 1:2:4:4 notitle palette