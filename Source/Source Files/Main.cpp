using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "Class_Definitions.h"  //Header File containing class definitions

void PrintToFile(Surface* surface)
{
	std::ofstream out("AFM.dat");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to AFM.dat

	surface->Print();

	std::cout.rdbuf(coutbuf); //reset to standard output again
}

int	main()
{
	int a1_cells = 10;
	int	a2_cells = 10;
	int a3_cells = 1;
	double a1[3] = { 3.0 / 2.0  , 3.0 / (2.0*sqrt(3.0)),  0.0 };
	double a2[3] = { 3.0 / 2.0  , -3.0 / (2.0*sqrt(3.0)),  0.0 };
	double a3[3] = { 0.25, 0.25, 1 };
	
				//x    y      z    r    m   eps sigma;
	Atom Carbon1(0.00, 0.00, 0.0, 1.0, 1.0,	0.0, 0.0);	 
	Atom Carbon2(0.33, 0.33, 0.0, 1.0, 1.0,	0.0, 0.0);	


	Unit_Cell* simple_cell = new Unit_Cell();
		simple_cell->Add_Atom(Carbon1);
		simple_cell->Add_Atom(Carbon2);

	Lattice* Square = new Lattice(a1, a2, a3);

	Surface* surface = new Surface(Square, simple_cell, a1_cells, a2_cells, a3_cells);

	PrintToFile(surface);

	return(0);
}
