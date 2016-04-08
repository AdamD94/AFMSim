#include "Var.h"


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
	MoveFile(Tip_Filename, Defect_Filename, "Topology.dat", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.dat", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Surface.dat", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Locations.dat", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "VestaObj.dat", Setpoint, z_0, Ang, Defect_Present);

	MoveFile(Tip_Filename, Defect_Filename, "Surface.png", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Force_Curve.png", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Topology.png", Setpoint, z_0, Ang, Defect_Present);
	MoveFile(Tip_Filename, Defect_Filename, "Derivative.png", Setpoint, z_0, Ang, Defect_Present);

	MoveFile(Tip_Filename, Defect_Filename, "Topology.xyz", Setpoint, z_0, Ang, Defect_Present);

}

void Plot()
{
	Gnuplot gp; 
	gp << "load 'Force_Curve.plt' \n";
	gp << "load 'Surface.plt' \n";
	gp << "load 'Topology.plt' \n";
	gp << "load 'Derivative.plt' \n";
}




int	main(int argc, char *argv[])
{

	if (argc < 2)
	{
		cout << "Incorrect number of arguments, call AFMSim <Tip Filename> <Defect Filename (optional)>" << endl;
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

	Atom* Carbon1 = new Atom(0.0, 0.0, 0.0, "C");		// Atoms to be added to the basis, coordinates are in terms of a1, a2, a3
	Atom* Carbon2 = new Atom((2.0/3.0), (2.0/3.0), 0.0, "C");
	
	Tip->ResetOrigin();
	Tip->Move(0 - Tip->r[0], 0 - Tip->r[1], z_0 - Tip->r[2]);
	
	Unit_Cell* simple_cell = new Unit_Cell();	// Unit cell to be tiled over the crystal
	simple_cell->Add_Atom(Carbon1);			// Constiuent atoms of unit cells
	simple_cell->Add_Atom(Carbon2);

	Surface* surface = new Surface(a1, a2, a3, a1_cells, a2_cells, a3_cells, simple_cell);	//Surface generation, unit cell is tiled over space in the directions defined by a1,a2,a3
	if (argc > 2)
	{
		cout << "Enter defect rotation:" << endl;
		cin>>Ang;
		Defect_Filename = argv[2];
		Defect->Import(Defect_Filename);
		Defect->Move(((a1[0] * a1_cells + a2[0]*a2_cells)/2)-Tip->r[0], ((a1[1] * a1_cells + a2[1] * a2_cells) / 2) - Tip->r[1], 1 - Tip->r[0]);	
		Defect->RotateAboutZ(Ang);
		Defect_Present = 1;
		surface->Add_Defect(Defect);
	}

	surface->Print();

	cout << "Computing" << endl;
	cout << "Cell count: " << surface->cell_count << endl;

	surface->Print();
	Tip->Print_Atoms();

/*
	cout<<"Calculating Frequency Shift"<<endl;
	Cantilever* Cant = new Cantilever(400E-6, 30E-6, 3.7E-6, 200 * pow(10, 9), 3184, Tip, surface);
	Cant->FM_Scan(surface, 1);
	gp << "load 'Transform.plt' \n";
	
	cout<<"Frequency calculation complete"<<endl;
*/

	surface->TipHeight(Tip, Setpoint, ZStep, FRes, 5);
	cout << "Topology complete" << endl;

	write_gsf("Topology.dat");

	surface->ForceCurve(Tip, z_0, z_0 + 1);
	cout << "Force curve complete" << endl;


	surface->SurfaceForce(Tip, z_0, 5);
	cout << "Surface force complete" << endl;


	Plot();

	MoveFiles();

	clock_t toc = clock(); // Read in current clock cycle, subtract and divide by clock frequency for time elapsed in seconds
	cout << "Simulation Complete \nElapsed: " << (double)(toc - tic) / CLOCKS_PER_SEC << "  seconds" << endl;
	return(0);
}