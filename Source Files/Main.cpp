using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <direct.h>
#include "Class_Definitions.h"  //Header File containing class definitions

void MoveFiles(string filename, double setpoint, double z_0, bool defect)
{
	string Directory = filename.substr(0, filename.find_last_of("\\"));
	string Folder = filename.substr(0, filename.find_last_of("."));
	Folder.append("_");
	Folder.append(to_string(setpoint).substr(0,5));
	Folder.append("_");
	Folder.append(to_string(z_0).substr(0, 5));
	if (defect == 1)
		Folder.append("_Defect");

	_mkdir(Folder.c_str());

	string OldFile = Directory;
	OldFile.append("\\Topology.png");
	string NewFile = Folder;
	NewFile.append("\\Topology.png");
	rename(OldFile.c_str(), NewFile.c_str());

	OldFile = Directory;
	OldFile.append("\\Topology.dat");
	NewFile = Folder;
	NewFile.append("\\Topology.dat");
	rename(OldFile.c_str(), NewFile.c_str());

	OldFile = Directory;
	OldFile.append("\\Force_Curve.png");
	NewFile = Folder;
	NewFile.append("\\Force_Curve.png");
	rename(OldFile.c_str(), NewFile.c_str());

	OldFile = Directory;
	OldFile.append("\\Force_Curve.dat");
	NewFile = Folder;
	NewFile.append("\\Force_Curve.dat");
	rename(OldFile.c_str(), NewFile.c_str());

	OldFile = Directory;
	OldFile.append("\\Surface.png");
	NewFile = Folder;
	NewFile.append("\\Surface.png");
	rename(OldFile.c_str(), NewFile.c_str());

	OldFile = Directory;
	OldFile.append("\\Surface.dat");
	NewFile = Folder;
	NewFile.append("\\Surface.dat");
	rename(OldFile.c_str(), NewFile.c_str());

	OldFile = Directory;
	OldFile.append("\\Locations.dat");
	NewFile = Folder;
	NewFile.append("\\Locations.dat");
	rename(OldFile.c_str(), NewFile.c_str());;

	OldFile = Directory;
	OldFile.append("\\Tip.dat");
	NewFile = Folder;
	NewFile.append("\\Tip.dat");
	rename(OldFile.c_str(), NewFile.c_str());;
}

int	main(int argc, char *argv[])
{

	double a = 1.4;
	double c = 6.71;
	double z_0 = 4;
	double ZStep = 1;
	double FRes = 0.0005; //Tolerance on force measurement
	double Setpoint = 0; //Force to be followed by tip nanoNewtons

	int a1_cells = 20;		// Number of cells to be generated in the a1, a2 and a3 directions
	int	a2_cells = 20;		
	int a3_cells = 1;	

	double a1[3] = {a* 3.0 / 2.0	, a * sqrt(3.0) / 2.0	,a*  0.0 };	// Lattice vectors as in crystallagraphy
	double a2[3] = {a* 3.0 / 2.0	, a *-sqrt(3.0) / 2.0	,a*  0.0 };	
	double a3[3] = {a				, 0						,c* -1.0 };	
	bool Defect_Present = 0;
	string Defect_Filename;
	string Tip_Filename;
	VestaObject* Defect= new VestaObject(0,0,0);

	if (argc < 2)
	{
		cout << "Incorrect number of arguments, drop tip onto .exe" << endl;
		cin.ignore();
		return(0);
	}


	Tip_Filename = argv[1];
	cout << "Tip: " << Tip_Filename << endl;

	cout << "Enter Force Setpoint (nN):" << endl; //Height of tip above sample is set by user
	cin >> Setpoint;

	cout << "Enter Height for surface force (Angstroms):" << endl; //Height of tip above sample is set by user
	cin >> z_0;

	clock_t tic = clock();	// Read in current clock cycle

	VestaObject* Tip = new VestaObject(0, 0, z_0);
	Tip->Import(Tip_Filename);
	Atom* Carbon1 = new Atom(0.0, 0.0, 0.0);		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom* Carbon2 = new Atom((2.0/3.0), (2.0/3.0), 0.0);
	
	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
	simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
	simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3
	
	if (argc > 2)
	{
		cout <<" defect "<< endl;
		Defect_Filename = argv[2];
		Defect->Import(Defect_Filename);
		Defect->Move(40.0, 0.0, 1.0);
		Defect_Present = 1;
		surface->Add_Defect(Defect);
	}

	cout << "Computing" << endl;
	cout << "Cell count: " << surface->cell_count << endl;

	surface->Print();
	Tip->Print_Atoms();

	surface->TipHeight(Tip, Setpoint, ZStep, FRes, 3);
	cout << "Topology complete" << endl;
	system("Topology.plt");
	
	surface->ForceCurve(Tip, z_0, z_0+1);
	cout << "Force curve complete" << endl;
	system("XZU.plt");

	surface->SurfaceForce(Tip, 1);
	cout << "Surface force complete" << endl;
	system("XYU.plt");

	MoveFiles(Tip_Filename, Setpoint, z_0, Defect_Present);

	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	cin.ignore();
	cin.ignore(); //Wait for input to terminate execution

	return(0);
}



