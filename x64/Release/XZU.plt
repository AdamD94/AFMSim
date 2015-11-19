#xy plot with potential 

set term pngcairo enhanced color crop size 1024,798 
set output 'Force Curve.png'

set palette model RGB
set palette model RGB defined (0 "dark-blue", 2 "green", 3 "yellow", 4 "red", 5 "black" )

unset size
unset origin

set xlabel 'x [Å]'
set ylabel 'z [Å]'
set cblabel 'F [nN]'

plot 'Force_Curve.dat' using 1:3:4 notitle palette