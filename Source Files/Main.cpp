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

	cout << "x "<<"y "<<"z "<<"r "<<"m "<< endl;

	surface->Print();

	std::cout.rdbuf(coutbuf); //reset to standard output again
}

void LJToFile(Surface* surface, Atom* Tip, int a1_cells ,int a2_cells, double a1[], double a2[])
{
	int res = 65;
	std::ofstream out("AFM.dat");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to AFM.dat

	cout << "x " << "y " << "U "<< endl;

	for (int j = (-res); j < a2_cells*res; j++)
	{
		for (int i = (-res); i < a1_cells*res; i++)
		{
			Tip->r[0] = (i*a1[0] + j*a2[0])/ res;
			Tip->r[1] = (i*a1[1] + j*a2[1])/ res;
			//Tip->XYPrint(surface->LJPot(Tip));
			Tip->XYPrint(surface->PerpLJForce(Tip));
		}
	}
	cout << endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
}


int	main()
{
	clock_t tic = clock();	// Read in current clock cycle
	int a1_cells = 10;		// Number of cells to be generated in the a1 direction
	int	a2_cells = 10;		// Number of cells to be generated in the a2 direction
	int a3_cells = 2;		// Number of cells to be generated in the a3 direction
	double a = 0.14;

	double a1[3] = {a* 3.0 / 2.0	, a *sqrt(3.0) / 2	,a*  0.0 };	// Lattice vector a1 as in crystallagraphy
	double a2[3] = {a* 3.0 / 2.0	,-a *sqrt(3.0) / 2	,a*  0.0 };	// Lattice vector a1 as in crystallagraphy
	double a3[3] = {a* 0.5			, a *sqrt(3.0)		,a*  -1.0 };	// Lattice vector a1 as in crystallagraphy

				//a1   a2    a3    r    m   element;
	Atom Carbon1(0  , 0  , 0, 0.4, 12, "C12");		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom Carbon2(1/3, 1/3, 0, 0.4, 12, "C12");	


	Atom* Tip=new Atom(0, 0, 5, 0.4, 12, "C12");	// Atom for the tip, coordinates are absolute

	Lattice* Square = new Lattice(a1, a2, a3);	// Lattice used to define shape of unit cells

	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
		simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
		simple_cell->Add_Atom(Carbon2);


	Surface* surface = new Surface(Square, simple_cell, a1_cells, a2_cells, a3_cells);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3

	LJToFile(surface, Tip, a1_cells, a2_cells, a1, a2);		//Print the LJ-potential to a file over extent of crystal

	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	
	cin.ignore(); //Wait for input to terminate execution

	return(0);
}
