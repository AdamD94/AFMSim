#xy plot with potential 

set term pngcairo size 1024,798
set output 'GNU_OUT.png'

set palette rgbformulae 22,13,-31
unset size
unset origin
set size ratio -1  
set xlabel 'x [nm]'
set ylabel 'y [nm]'
set cblabel 'F [nN]'

plot 'AFM.dat' using 1:2:4 notitle palette