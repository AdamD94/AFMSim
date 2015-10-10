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
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

	surface->Print_Lattice_Atom();

	std::cout.rdbuf(coutbuf); //reset to standard output again
}

int	main()
{
	int a1_cells = 10;
	int	a2_cells = 10;
	int a3_cells = 1;
	double a1[3];  
	double a2[3]; 
	double a3[3]; 
	double offset[3] = {0.5, 0.5, 0.5};
	
	a1[0] = 3.0 / 2.0;
	a1[1] = (3.0 / (2.0*sqrt(3.0)));
	a1[2] = 0.0;

	a2[0] = 3.0 / 2.0;
	a2[1] = (-3.0 / (2.0*sqrt(3.0)));
	a2[2] = 0.0;

	a3[0] = 0.0;
	a3[1] = 0.0;
	a3[2] = 1.0;

	Atom Carbon1(0.0, 0.0, 0.0, 1.0, 1.0, 0, 0);  //x,y,z,r,m, epsilon, sigma;
	Atom Carbon2(0.33, 0.33, 0.0, 1.0, 1.0, 0, 0);  //x,y,z,r,m, epsilon, sigma;


	Unit_Cell* simple_cell = new Unit_Cell();

	simple_cell->Add_Atom(Carbon1);
	simple_cell->Add_Atom(Carbon2);

	Lattice* Square = new Lattice(a1, a2, a3);
//	Square->next = new Lattice(a1, a2, a3, offset);

	Surface* surface = new Surface(Square, simple_cell, a1_cells, a2_cells, a3_cells);

	PrintToFile(surface);



	return(0);
}
