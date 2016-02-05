using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include "Class_Definitions.h"  //Header File containing class definitions


int	main(int argc, char* argv[])
{
	int a1_cells = 14;		// Number of cells to be generated in the a1, a2 and a3 directions
	int	a2_cells = 14;		
	int a3_cells = 1;	

	double a = 1.4;
	double c = 6.71;
	double aPt = 2.27;
	double z_0 = 4;

	double a1[3] = {a* 3.0 / 2.0	, a * sqrt(3.0) / 2.0	,a*  0.0 };	// Lattice vectors as in crystallagraphy
	double a2[3] = {a* 3.0 / 2.0	, a *-sqrt(3.0) / 2.0	,a*  0.0 };	
	double a3[3] = {a				, 0						,c* -1.0 };	

	double ZRes = 0.1; // vertical resolution in angstroms
	double FRes = 0.0005; //Tolerance on force measurement
	double Setpoint = 0; //Force to be followed by tip nanoNewtons

	
	if (argc != 2) //Error handling
	{
		cout << "ERROR: Wrong amount of arguments!" << endl;
		cout << "\nDrag .vesta file onto executable \nwith .xyz file with same name in same directory" << endl;
		cin.ignore();
		return 0;
	}

	string filename = argv[1];

	cout << "Enter Force Setpoint (nN):" << endl; //Height of tip above sample is set by user
	cin >> Setpoint;

	clock_t tic = clock();	// Read in current clock cycle

	Tip* Tip1 = new Tip(0, 0, z_0);
		Tip1->ImportTip(filename);
		Tip1->Print_Atoms();
				
	Atom* Carbon1 = new Atom(0.0, 0.0, 0.0);		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom* Carbon2 = new Atom((2.0/3.0), (2.0/3.0), 0.0);

	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
		simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
		simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3

	cout << "Computing" << endl;

	cout << "Cell count: "<<surface->cell_count << endl;
	surface->Print();
	
	surface->TipHeight(Tip1, Setpoint, ZRes, FRes);

/*
	surface->ForceCurve(Tip1, z_0, z_0+1);
	cout << "Force curve complete" << endl;

	surface->SurfaceForce(Tip1);
	cout << "Surface force complete" << endl;
*/
	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	cin.ignore();
	cin.ignore(); //Wait for input to terminate execution

	return(0);
}
