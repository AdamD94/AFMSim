
set term pngcairo size 1024,798
set output 'GNU_OUT.png'
set palette rgbformulae 22,13,-31
set xlabel  'x [nm]'
set ylabel  'z [nm]'
set cblabel 'F [nN]'
plot 'Force_Curve.dat' using 1:3:4 notitle palette