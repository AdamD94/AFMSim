#include "Var.h"
#include "AFMSim.h"

//Generate folder strucure, rename and store generated data
void MoveFile(string Tipname, string Defectname, string Filename, double ID1, double ID2, double ID3, bool DefectPresent)
{
	Defectname=Defectname.substr(Defectname.find_last_of("/")+1,1000);
	Tipname=Tipname.substr(Tipname.find_last_of("/")+1,1000);

	Defectname=Defectname.substr(0,Defectname.find_last_of("."));
	Tipname=Tipname.substr(0,Tipname.find_last_of("."));

	string Directory = "Results";
	create_directory(Directory);
	Directory.append("/");

	if(DefectPresent)
	{
		Directory.append(Defectname);
		create_directory(Directory);
		Directory.append("/");
	}

	Directory.append(Tipname);
	Directory.append("/");
	create_directory(Directory);

	string Newfile = Directory;

	Newfile.append(to_string(ID1).substr(0,5));
	Newfile.append("_");
	Newfile.append(to_string(ID2).substr(0,4));
	Newfile.append("_");
	Newfile.append(to_string(ID3).substr(0,4));
	Newfile.append("_");
	Newfile.append(Filename);

	boost::filesystem::rename(Filename,Newfile);
}

void MoveFiles()
{
	MoveFile(Tip_Filename, Defect_Filename, "Topology.dat",		Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.dat", 	Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Surface.dat", 		Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Locations.dat", 	Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "VestaObj.dat", 	Setpoint, z_0, Ang, Defect_Present);

	MoveFile(Tip_Filename, Defect_Filename, "Surface.png", 		Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.png", 	Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Combination.png",	Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Topology.png", 	Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Derivative.png",	Setpoint, z_0, Ang, Defect_Present);


	MoveFile(Tip_Filename, Defect_Filename, "Topology.xyz", 	Setpoint, z_0, Ang, Defect_Present);
}

void Plot()
{
	Gnuplot gp;
	gp << "load 'Combination.plt' \n";
	gp << "load 'Force_Curve.plt' \n";
	gp << "load 'Surface.plt' \n";
	gp << "load 'Topology.plt' \n";
	gp << "load 'Derivative.plt' \n";
}



int	main(int argc, char *argv[])
{
	if (argc < 2) //Check for simulation call without parmeters
	{
		cout << "Error: Incorrect number of arguments\nCall AFMSim <TipFilename> <DefectFilename>(optional)" << endl;
		return(0);
	}

	Tip_Filename = argv[1];
	cout << "Tip: " << Tip_Filename << endl;

	cout << "Enter Force Setpoint (nN):" << endl; 					// Force to be followed by tip to generate topology
	cin >> Setpoint;

	cout << "Enter Height for surface force (Angstroms):" << endl; 	// Initial height of tip apex above sample is set by user
	cin >> z_0;

	VestaObject* Tip = new VestaObject(0, 0, z_0);
	Tip->Import(Tip_Filename);

	Tip->ResetOrigin();
	Tip->Move(0 - Tip->r[0], 0 - Tip->r[1], z_0 - Tip->r[2]);

	cout 	<< "Functionalise tip with oxygen atom? [y/n]" << endl;
	cin 	>> Functionalise;

	if(strcmp(Functionalise.c_str(),"Y")==0||strcmp(Functionalise.c_str(),"y")==0)	// Add oxygen atom 1 angstrom below tip apex
	{																				// if functionalise option is selected
		
		Atom* Oxygen = new Atom(0,0,-1,"O");
		Tip->Add_Atom(Oxygen);
		Tip->ResetOrigin();
		Tip->Move(0 - Tip->r[0], 0 - Tip->r[1], z_0 - Tip->r[2]);
	}

	// Atoms to be added to the basis, coordinates are fractional in terms of a1, a2, a3
	Atom* Carbon1 = new Atom(0.0, 		0.0, 			0.0, "C");		
	Atom* Carbon2 = new Atom((2.0/3.0), (2.0/3.0), 	0.0, "C");

	Unit_Cell* simple_cell = new Unit_Cell();		// Unit cell to be tiled over the crystal
	simple_cell->Add_Atom(Carbon1);					// Constiuent atoms of unit cells
	simple_cell->Add_Atom(Carbon2);

	// Surface generation: unit cell is tiled over space in the directions defined by a1,a2,a3
	Surface* 		surface	=	new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell); 	
	VestaObject* 	Defect 	= 	new VestaObject(0, 0, 0);														

	if (argc > 2)
	{
		cout << "Enter defect rotation:" << endl;
		cin>>Ang;
		Defect_Filename = argv[2];
		Defect->	Import(Defect_Filename);
		Defect->	Move(((a1[0] * a1_cells + a2[0]*a2_cells)/2)-Tip->r[0], ((a1[1] * a1_cells + a2[1] * a2_cells) / 2) - Tip->r[1], 1 - Tip->r[0]);
		Defect->	RotateAboutZ(Ang);
		surface->	Add_Defect(Defect);
		Defect_Present = 1;
	}

	clock_t tic = clock();	// Read in current clock cycle
	cout << "Computing" << endl;
	cout << "Cell count: " << surface->cell_count << endl;

	surface->Print();
	Tip->Print_Atoms();

	surface->TipHeight(Tip, Setpoint, ZStep, FRes, 5); 	// generates AFM-like topology
	cout << "Topology complete" << endl;

	write_xyz("Topology.dat");							// Writes gwyddion-compatible file

	surface->ForceCurve(Tip, z_0, z_0 + 1, 1);				// Generate x-z force data
	cout << "Force curve complete" << endl;

	surface->SurfaceForce(Tip, z_0, 1);					// Generate x-y force data
	cout << "Surface force complete" << endl;

	Plot();						// Call Gnuplot

	MoveFiles();				// Store results in folder structure & generate folder structure if necessary

	clock_t toc = clock();		// Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	return(0);
}