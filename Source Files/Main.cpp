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
	clock_t tic = clock();	// Read in current clock cycle

	int a1_cells = 12;		// Number of cells to be generated in the a1, a2 and a3 directions
	int	a2_cells = 12;		
	int a3_cells = 1;	

	double a = 1.4;
	double c = 6.71;
	double aPt = 2.27;

	double a1[3] = {a* 3.0 / 2.0	, a * sqrt(3.0) / 2.0	,a*  0.0 };	// Lattice vectors as in crystallagraphy
	double a2[3] = {a* 3.0 / 2.0	, a *-sqrt(3.0) / 2.0	,a*  0.0 };	
	double a3[3] = {a			, 0						,c* -1.0 };	

	double e = 3.534*pow(10, -21);	//epsilon	[Joules]	(LJ Parameter) 
	double s = 2.95;				//sigma		[Angstroms] (LJ Parameter) 

	Atom* atom1 = new Atom( 0,   0,   0,   e, s);	// Atom for the tip, coordinates are absolute
	Atom* atom2 = new Atom(aPt / 4,   aPt/4,   aPt/5, e, s);
	Atom* atom3 = new Atom( aPt, 0,   aPt, e, s);
	Atom* atom4 = new Atom(-aPt, 0,   aPt, e, s);
	Atom* atom5 = new Atom( 0,   aPt, aPt, e, s);
	Atom* atom6 = new Atom( 0,  -aPt, aPt, e, s);

		Tip* Tip1 = new Tip(0, 0, 3);
			Tip1->Add_Atom(atom1);
			Tip1->Add_Atom(atom2);
			Tip1->Add_Atom(atom3);
			Tip1->Add_Atom(atom4);
			Tip1->Add_Atom(atom5);
			Tip1->Add_Atom(atom6);
	
				
	Atom* Carbon1 = new Atom(0.0, 0.0, 0.0);		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom* Carbon2 = new Atom((2.0/3.0), (2.0/3.0), 0.0);

	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
		simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
		simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3
	
	surface->ForceCurve(Tip1);
	surface->SurfaceForce(Tip1);
	surface->Print();

	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	cin.ignore(); //Wait for input to terminate execution

	return(0);
}
