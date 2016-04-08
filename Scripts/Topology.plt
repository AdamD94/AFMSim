#Tip Height Plot
unset xrange
unset yrange
unset zrange
unset cbrange

set term pngcairo enhanced color size 1024,798 crop font 'Veranda, 18'
set output 'Topology.png'
  
set palette model RGB
set palette rgb 34,35,36; 
set palette functions gray, gray, gray


set size ratio -1
set origin -0.033, 0

unset cblabel

set cblabel 'z [Å]'

set lmargin 0
set tmargin 0
set bmargin 0

unset xtics
unset ytics
unset border

set arrow 3 from 48,-8 to 53,-8 nohead size 1, 90 lw 3 lt -1 front 
set label "5Å" at 50,-8.75

set pm3d map 
splot 'Topology.dat' using 1:2:3:3 notitle palette