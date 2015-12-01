#xy plot with potential 

set term png 
set output 'Surface.png'
  
set palette model RGB
set palette model RGB defined (	0 "dark-blue",2 "green",3 "yellow",4 "red",5 "black")

set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'F [nN]'

set size ratio -1

set pm3d map interpolate 2,2
splot 'Surface.dat' using 1:2:4 notitle