#include "Main.h"



void MoveFile(string Tipname, string Defectname, string Filename, double ID1, double ID2, bool DefectPresent)
{
	std::string Directory = Tipname.substr(0, Tipname.find_last_of("\\"));
	std::string Folder = Tipname.substr(0, Tipname.find_last_of("."));
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

	create_directory(Folder.c_str());

	std::string OldFile = Directory;
	OldFile.append("\\");
	OldFile.append(Filename);
	std::string NewFile = Folder;
	NewFile.append("\\");
	NewFile.append(Filename);
	rename(OldFile.c_str(), NewFile.c_str());
}



void MoveFiles()
{
	MoveFile(Tip_Filename, Defect_Filename, "Topology.png", Setpoint, z_0, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Derivative.png", Setpoint, z_0, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Topology.dat", Setpoint, z_0, Defect_Present);

	MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.png", Setpoint, z_0, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.dat", Setpoint, z_0, Defect_Present);

	MoveFile(Tip_Filename, Defect_Filename, "Surface.png", Setpoint, z_0, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Surface.dat", Setpoint, z_0, Defect_Present);

	MoveFile(Tip_Filename, Defect_Filename, "Locations.dat", Setpoint, z_0, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "VestaObj.dat", Setpoint, z_0, Defect_Present);
}


int	main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "Incorrect number of arguments, drop tip onto .exe" << endl;
		std::cin.ignore();
		return(0);
	}
	
	Tip_Filename = argv[1];
	std::cout << "Tip: " << Tip_Filename << endl;

	std::cout << "Enter Force Setpoint (nN):" << endl; //Height of tip above sample is set by user
	std::cin >> Setpoint;

	std::cout << "Enter Height for surface force (Angstroms):" << endl; //Height of tip above sample is set by user
	std::cin >> z_0;

	clock_t tic = clock();	// Read in current clock cycle

	VestaObject* Tip = new VestaObject(0, 0, z_0);
	Tip->Import(Tip_Filename);

	Atom* Carbon1 = new Atom(0.0, 0.0, 0.0);		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom* Carbon2 = new Atom((2.0/3.0), (2.0/3.0), 0.0);
	Atom* TipCarbon = new Atom(0.0, 0.0, -3.0, 1.0456*pow(10, -21), 4.9);		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	
	Tip->ResetOrigin();
	Tip->Add_Atom(TipCarbon);
	Tip->ResetOrigin();
	Tip->Move(0 - Tip->r[0], 0 - Tip->r[1], z_0 - Tip->r[2]);
	
	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
	simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
	simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3
	if (argc > 2)
	{
		std::cout << " defect " << endl;
		Defect_Filename = argv[2];
		Defect->Import(Defect_Filename);
		Defect->Move(((a1[0] * a1_cells + a2[0]*a2_cells)/2)-Tip->r[0], ((a1[1] * a1_cells + a2[1] * a2_cells) / 2) - Tip->r[1], 1 - Tip->r[0]);
		Defect_Present = 1;
		surface->Add_Defect(Defect);
	}

	surface->Print();

	std::cout << "Computing" << endl;
	std::cout << "Cell count: " << surface->cell_count << endl;

	Cantilever* Cant = new Cantilever(400E-6, 30E-6, 3.7E-6, 200 * pow(10, 9), 3184, Tip, surface);
	Cant->FM_Scan(surface, 1);
	
	surface->Print();
	Tip->Print_Atoms();

	surface->TipHeight(Tip, Setpoint, ZStep, FRes, 3);
	cout << "Topology complete" << endl;

	surface->ForceCurve(Tip, z_0, z_0 + 1);
	cout << "Force curve complete" << endl;
		
	surface->SurfaceForce(Tip, z_0, 3);
	cout << "Surface force complete" << endl;
		
	MoveFiles();
	

	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	std::cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	std::cin.ignore();
	std::cin.ignore(); //Wait for input to terminate execution

	return(0);
}



