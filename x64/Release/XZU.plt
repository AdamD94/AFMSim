
set term pngcairo size 1024,798
set output 'Force_Curve.png'
set palette rgbformulae 22,13,-31
set xlabel  'x [Å]'
set ylabel  'z [Å]'
set cblabel 'F [nN]'
plot 'Force_Curve.dat' using 1:3:4 notitle palette