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

void LJToFile(Surface* surface, Atom* atom_in, int a1_cells ,int a2_cells, double a1[], double a2[])
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
			atom_in->r[0] = (i*a1[0] + j*a2[0])/ res;
			atom_in->r[1] = (i*a1[1] + j*a2[1])/ res;
			atom_in->XYPrint(surface->LJPot(atom_in));
		}
	}
	cout << endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
}

void LJToFile(Surface* surface, Tip* Tip_in, int a1_cells, int a2_cells, double a1[], double a2[])
{
	int res = 65;
	int k = 0;
	double Average = 0;
	std::ofstream out("AFM.dat");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to AFM.dat

	cout << "x " << "y " << "z " << "U " << endl;

	for (int j = (((a2_cells/2)-1)*res); j < (((a2_cells / 2) + 1)*res); j++)
	{
		for (int i = (((a1_cells / 2) - 1)*res); i < (((a1_cells / 2) + 1)*res); i++)
		{
			Tip_in->MoveTip((i*a1[0] + j*a2[0]) / res - Tip_in->r[0], (i*a1[1] + j*a2[1]) / res - Tip_in->r[1], 0);
			Average+=surface->LJForce(Tip_in);
			k++;
		}
	}

	Average = Average / k;

	for (int j = (-res); j < a2_cells*res; j++)
	{
		for (int i = (-res); i < a1_cells*res; i++)
		{
			Tip_in->MoveTip((i*a1[0] + j*a2[0]) / res - Tip_in->r[0], (i*a1[1] + j*a2[1]) / res - Tip_in->r[1], 0);
			Tip_in->Print(surface->LJForce(Tip_in)-Average);
		}
	}
	cout << endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
}

void ForceCurve(Surface* surface, Atom* Atom_in, int a1_cells, int a2_cells, double a1[], double a2[])
{
	double zstep = 0.0005;
	double xstep = 0.001;
	double Average = 0;
	double XMax = (((a1_cells / 2) + 0.5)*a1[0] + ((a2_cells / 2) + 0.5)*a2[0]);
	double XMin = (((a1_cells / 2) - 1.5)*a1[0] + ((a2_cells / 2) - 1.5)*a2[0]);
	double YMax = abs(a1[1]) + abs(a2[1]);
	double YMin = 0;
	double ZMax = 0.4;
	double ZMin = 0.3;
	int i=0;

	std::ofstream out("Force_Curve.dat");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to AFM.dat

	cout << "x " << "y " << "U " << endl;

	Atom_in->r[0] = XMin;
	Atom_in->r[1] = YMin;
	Atom_in->r[2] = ZMin;

	while (Atom_in->r[2] < ZMax)
	{
		Atom_in->r[1] = YMin;
		while (Atom_in->r[1] < YMax)
		{
			Atom_in->r[0] = XMin;
			while (Atom_in->r[0] < XMax)
			{
				Average += surface->LJForce(Atom_in);
				Atom_in->r[0] += xstep*10;
				i++;
			}
			Atom_in->r[1] += xstep*10;
		}
		Average = Average/i;

		Atom_in->r[0] = XMin;
		Atom_in->r[1] = YMin;

		while (Atom_in->r[0] < XMax)
		{
			Atom_in->Print(surface->LJForce(Atom_in)-Average);
			Atom_in->r[0] += xstep;
		}

		Atom_in->r[2] += zstep;
		Average = 0;
		i = 0;
		Atom_in->r[0] = XMin;
	}

	cout << endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
}

void ForceCurve(Surface* surface, Tip* Tip_in, int a1_cells, int a2_cells, double a1[], double a2[])
{
	double zstep = 0.0005;
	double xstep = 0.001;
	double Average = 0;
	double XMax = (((a1_cells / 2) + 0.5)*a1[0] + ((a2_cells / 2) + 0.5)*a2[0]);
	double XMin = (((a1_cells / 2) - 1.5)*a1[0] + ((a2_cells / 2) - 1.5)*a2[0]);
	double YMax = abs(a1[1]) + abs(a2[1]);
	double YMin = 0;
	double ZMax = 0.4;
	double ZMin = 0.3;
	int i = 0;

	std::ofstream out("Force_Curve.dat");
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to AFM.dat

	cout << "x " << "y " << "U " << endl;

	Tip_in->MoveTip(XMin - Tip_in->r[0], YMin - Tip_in->r[1], ZMin - Tip_in->r[2]);

	while (Tip_in->r[2] < ZMax)
	{
		Tip_in->MoveTip(0, YMin - Tip_in->r[1], 0);
		while (Tip_in->r[1] < YMax)
		{
			Tip_in->MoveTip(XMin - Tip_in->r[0], 0, 0);
			while (Tip_in->r[0] < XMax)
			{
				Average += surface->LJForce(Tip_in);
				Tip_in->MoveTip(xstep * 10, 0, 0);
				i++;
			}
			Tip_in->MoveTip(0, xstep * 10, 0);
		}
		Average = Average / i;

		Tip_in->MoveTip(XMin - Tip_in->r[0], YMin - Tip_in->r[1], 0);
		while (Tip_in->r[0] < XMax)
		{
			
			Tip_in->Print(surface->LJForce(Tip_in) - Average);
			Tip_in->MoveTip(xstep, 0, 0);
		}

		Tip_in->MoveTip(0, 0, zstep);
		Average = 0;
		i = 0;
		Tip_in->MoveTip(XMin - Tip_in->r[0], 0, 0);
	}
	cout << endl;
	std::cout.rdbuf(coutbuf); //reset to standard output again
}

int	main()
{
	clock_t tic = clock();	// Read in current clock cycle

	int a1_cells = 10;		// Number of cells to be generated in the a1, a2 and a3 directions
	int	a2_cells = 10;		
	int a3_cells = 2;	

	double a = 0.14;
	double a1[3] = {a* 3.0 / 2.0	, a * sqrt(3.0) / 2	,a*  0.0 };	// Lattice vectors as in crystallagraphy
	double a2[3] = {a* 3.0 / 2.0	, a *-sqrt(3.0) / 2	,a*  0.0 };	
	double a3[3] = {a* 0.5			, a * sqrt(3.0)		,a* -1.0 };	

				//a1   a2    a3    r    m   element
	Atom Carbon1(0  , 0  , 0, 0.4, 12, "C12");		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom Carbon2(1/3, 1/3, 0, 0.4, 12, "C12");	


	Atom* atom1 = new Atom( 0,  0, 0, 0.4, 12, "C12");	// Atom for the tip, coordinates are absolute
	Atom* atom2 = new Atom( 0,  0, a, 0.4, 12, "C12");
	Atom* atom3 = new Atom( a,  0, a, 0.4, 12, "C12");
	Atom* atom4 = new Atom(-a,  0, a, 0.4, 12, "C12");
	Atom* atom5 = new Atom( 0,  a, a, 0.4, 12, "C12");
	Atom* atom6 = new Atom( 0, -a, a, 0.4, 12, "C12");

	Tip* Tip1 = new Tip(0,0,0.3);
		Tip1->Add_Atom(atom1);
		Tip1->Add_Atom(atom2);
		Tip1->Add_Atom(atom3);
		Tip1->Add_Atom(atom4);
		Tip1->Add_Atom(atom5);
		Tip1->Add_Atom(atom6);

	Lattice* Square = new Lattice(a1, a2, a3);	// Lattice used to define shape of unit cells

	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
		simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
		simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(Square, simple_cell, a1_cells, a2_cells, a3_cells);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3

//	LJToFile(surface, Tip1, a1_cells, a2_cells, a1, a2);
	ForceCurve(surface, Tip1, a1_cells, a2_cells, a1, a2);	//Print the LJ-potential to a file over extent of crystal
	
	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	cin.ignore(); //Wait for input to terminate execution

	return(0);
}
