# AFMSim
This project is an atomic force microscopy simulation written using object-oriented C++.
It utilises a truncated many-body Lennard-Jones potential.

It was written as part of a final year physics project and is capable of producing images matching closely with experiment.
Forces can also be calculated in xy, xz, or yz planes.

AFM tips and surfaces can be imported from a crystallographic program, Vesta: http://jp-minerals.org/vesta/en/ 

Data is plotted with gnuplot automatically as long as gnuplot is installed: http://www.gnuplot.info/

Topology data is saved in a format that can be oppened with Gwyddion http://gwyddion.net/

Note, this program has dependencies on the boost library.
This library can be installed on linux with the following command.
sudo apt-get install libboost-all-dev
