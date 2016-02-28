#include "Main.h"

void MoveFile(string Tipname, string Defectname, string Filename, double ID1, double ID2, bool DefectPresent)
{

	string Directory = Tipname.substr(0, Tipname.find_last_of("\\"));
	string Folder = Tipname.substr(0, Tipname.find_last_of("."));
	Folder.append("_");
	if (DefectPresent)
	{
		Defectname = Defectname.substr(Defectname.find_last_of("\\") + 1, Defectname.find_last_of("."));
		Defectname = Defectname.substr(0, Defectname.find_last_of("."));
		Folder.append(Defectname);
		Folder.append("_");
	}
	Folder.append(to_string(ID1).substr(0,6));
	Folder.append("_");
	Folder.append(to_string(ID2).substr(0,6));

	_mkdir(Folder.c_str());

	string OldFile = Directory;
	OldFile.append("\\");
	OldFile.append(Filename);
	string NewFile = Folder;
	NewFile.append("\\");
	NewFile.append(Filename);
	rename(OldFile.c_str(), NewFile.c_str());
}

int	main(int argc, char *argv[])
{
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
	Atom* TipCarbon = new Atom(0.0, 0.0, -3.0, 1.0456*pow(10, -21), 4.9);		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	
	Tip->ResetOrigin();
	//Tip->Add_Atom(TipCarbon);
	Tip->ResetOrigin();
	Tip->Move(0 - Tip->r[0], 0 - Tip->r[1], z_0 - Tip->r[2]);
	
	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
	simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
	simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3
	if (argc > 2)
	{
		cout << " defect " << endl;
		Defect_Filename = argv[2];
		Defect->Import(Defect_Filename);
		Defect->Move(((a1[0] * a1_cells + a2[0]*a2_cells)/2)-Tip->r[0], ((a1[1] * a1_cells + a2[1] * a2_cells) / 2) - Tip->r[1], 1 - Tip->r[0]);
		Defect_Present = 1;
		surface->Add_Defect(Defect);
	}

	cout << "Computing" << endl;
	cout << "Cell count: " << surface->cell_count << endl;

	while (Setpoint <= 1)
	{
		surface->Print();
		Tip->Print_Atoms();

		surface->TipHeight(Tip, Setpoint, ZStep, FRes, 3);
		cout << "Topology complete" << endl;
		system("Topology.plt");
		system("Derivative.plt");
		
		surface->ForceCurve(Tip, z_0, z_0 + 1);
		cout << "Force curve complete" << endl;
		system("XZU.plt");
		
		surface->SurfaceForce(Tip, z_0, 3);
		cout << "Surface force complete" << endl;
		system("XYU.plt");
		
		MoveFile(Tip_Filename, Defect_Filename, "Topology.png"  ,	Setpoint, z_0, Defect_Present);
		MoveFile(Tip_Filename, Defect_Filename, "Derivative.png",	Setpoint, z_0, Defect_Present);
		MoveFile(Tip_Filename, Defect_Filename, "Topology.dat"  , Setpoint, z_0, Defect_Present);
		
		MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.png", Setpoint, z_0, Defect_Present);
		MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.dat", Setpoint, z_0, Defect_Present);
		
		MoveFile(Tip_Filename, Defect_Filename, "Surface.png", Setpoint, z_0, Defect_Present);
		MoveFile(Tip_Filename, Defect_Filename, "Surface.dat", Setpoint, z_0, Defect_Present);
		
		MoveFile(Tip_Filename, Defect_Filename, "Locations.dat", Setpoint, z_0, Defect_Present);
		MoveFile(Tip_Filename, Defect_Filename, "VestaObj.dat" , Setpoint, z_0, Defect_Present);
		Setpoint += 0.01;
		cout <<"Setpoint: "<< Setpoint << endl;
	}

	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	cin.ignore();
	cin.ignore(); //Wait for input to terminate execution

	return(0);
}



