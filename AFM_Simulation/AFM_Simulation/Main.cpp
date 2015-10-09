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
	int a3_cells = 5;
	double a1[3] = {1,0,0};
	double a2[3] = {0,1,0};
	double a3[3] = {0,0,1};
	double offset[3] = { 0.5, 0.0, 0.5 };
	
	Atom Carbon1(0, 0.0, 0.0, 1.0, 1.0);  //x,y,z,r,m;
	Atom Carbon2(0.5, 0.5, 0.0, 1.0, 1.0);  //x,y,z,r,m;


	Unit_Cell* simple_cell = new Unit_Cell();

	simple_cell->Add_Atom(Carbon1);
	simple_cell->Add_Atom(Carbon2);

	Lattice* Square = new Lattice(a1, a2, a3);
	Lattice* Square2 = new Lattice(a1, a2, a3, offset);
	//Square->next = Square2;

	Surface* surface = new Surface(Square, simple_cell, a1_cells, a2_cells, a3_cells);

	PrintToFile(surface);



	return(0);
}
