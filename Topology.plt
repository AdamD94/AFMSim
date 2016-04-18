#Tip Height Plot
set term pngcairo enhanced color size 1024,798 crop font 'Veranda, 18'
set output 'Topology.png'
  
set palette model RGB
set palette rgb 34,35,36; 
set palette functions gray, gray, gray

set size ratio -1
set origin -0.033, 0

unset colorbox
unset cblabel
set cblabel 'z [Å]'

unset xlabel
unset ylabel
unset xtics
unset ytics
unset border

set pm3d map interpolate 2,2
splot 'Topology.dat' using 1:2:3:3 notitle palette