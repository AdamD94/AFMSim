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
	int res = 100.0;
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
			Tip->XYPrint(surface->LJPot(Tip));	
		}
	}
	cout << endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
}


int	main()
{
	clock_t tic = clock();
	int a1_cells = 10;
	int	a2_cells = 10;
	int a3_cells = 1;
	double a1[3] = { 3.0 / 2.0  ,  3.0 / (2.0*sqrt(3.0)),  0.0 };
	double a2[3] = { 3.0 / 2.0  , -3.0 / (2.0*sqrt(3.0)),  0.0 };
	double a3[3] = { 0.5		,  3.0 / (2.0*sqrt(3.0)),  -1.0 };

				//x    y      z    r    m   element;
	Atom Carbon1(0.0 , 0.0 , 0.0, 0.4, 1.0, "C12");	 
	Atom Carbon2(0.33, 0.33, 0.0, 0.4, 1.0,	"C12");	

	Atom* Tip=new Atom(0.0, 0.0, 1.0, 1.0, 1.0, "C12");


	Unit_Cell* simple_cell = new Unit_Cell();
		simple_cell->Add_Atom(Carbon1);
		simple_cell->Add_Atom(Carbon2);

	Lattice* Square = new Lattice(a1, a2, a3);

	Surface* surface = new Surface(Square, simple_cell, a1_cells, a2_cells, a3_cells);

	LJToFile(surface, Tip, a1_cells, a2_cells, a1, a2);

	clock_t toc = clock();
	cout << "Simulation Complete \n Elapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	
	cin.ignore();

	return(0);
}
