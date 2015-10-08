using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "Class_Definitions.h"  //Header File containing class definitions


int	main()
{
	double x_max = 20; //20nm x extent
	double y_max = 20; //20nm y extent
	double a1[3] = { 1, 0, 0 };
	double a2[3] = { 0, 1, 0 };


	Atom Carbon(0.0, 0.0, 0.0, 1.0, 1.0);  //x,y,z,r,m;
	Unit_Cell* simple_cell = new Unit_Cell();
	simple_cell->Add_Atom(Carbon);
	Lattice* Square = new Lattice(a1, a2);
	Surface* surface = new Surface();

	surface->Generate_Surface(Square, simple_cell, x_max, y_max);

	std::ofstream out("AFM.dat");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	surface->Print();


	std::cout.rdbuf(coutbuf); //reset to standard output again

	return(0);
}
