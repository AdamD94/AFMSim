#xy plot with potential 

set term png 
set output 'Surface.png'
  
set palette model RGB
set palette model RGB defined (	0 "dark-blue",2 "green",3 "yellow",4 "red",5 "black")

unset size
unset origin
set size ratio -1

set xlabel 'x [Å]'
set ylabel 'y [Å]'
set cblabel 'F [nN]'

plot 'Surface.dat' using 1:2:4 notitle palette 