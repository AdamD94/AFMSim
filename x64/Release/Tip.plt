#xyz plot of probe tip 

  
unset cblabel
set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'z [Å]'

set view equal xyz

splot 'Tip.dat' using 1:2:3 notitle  pt 7 ps 15
