#ifndef Class_Definitions
#define Class_Definitions

class Atom
{
private:
	double temp;

public:
	double r[3];
	double	e;
	double	s;
	
	Atom* Next;

	Atom(double x_in, double y_in, double z_in)
	{
		r[0] = x_in;
		r[1] = y_in;
		r[2] = z_in;
		e = 0;
		s = 0;	
		Next = NULL;
	}

	Atom(double x_in, double y_in, double z_in, double e_in, double s_in)
	{
		r[0] = x_in;
		r[1] = y_in;
		r[2] = z_in;
		e = e_in;
		s = s_in;
		Next = NULL;
	}

	void Print()
	{
		cout << r[0] << " " << r[1] << " " << r[2] << "\n";
	}

	void Print(int var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << var << "\n";
	}

	void Print(double var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << var << "\n";
	}

	void Print(double var1, double var2)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << var1 << " " << var2 << "\n";
	}

	double Dist(Atom* Other_Atom)
	{
		return (double)  sqrt(pow((Other_Atom->r[2] - r[2]), 2) + pow((Other_Atom->r[1] - r[1]), 2) + pow((Other_Atom->r[0] - r[0]), 2));
	}

	double LJPot(Atom* Other_Atom)
	{	
		temp = Dist(Other_Atom);
		if (temp > (3 * s))
			return 0;

		const double LJ_Truncation = 4 * e*(pow((s / (3 * s)), 12) - pow((s / (3 * s)), 6));

		return (double)	 4 * e*(pow((s /temp), 12) - pow((s / temp), 6)) - LJ_Truncation;
	}

	double LJForce(Atom* Other_Atom)
	{
		temp = Dist(Other_Atom);
		if (temp > (3 * s))
			return 0;

		const double LJ_Truncation = (4.0 / (3*s))*e*(12 * pow((s / (3 * s)), 12) - 6 * pow((s / (3 * s)), 6))*pow(10, 19);

		return (double) (4.0 / temp)*e*(12 * pow((s / temp), 12) - 6 * pow((s / temp), 6))*pow(10, 19) - LJ_Truncation;
	}
	
	double PerpLJForce(Atom* Other_Atom)
	{
		return (double) LJForce(Other_Atom)*(abs(Other_Atom->r[2]-r[2]))/(Dist(Other_Atom));
	}

};

class VestaObject
{
private:

	double orientation[4][4];

public:

	Atom* atom;
	int atom_count=0;
	double r[3];


	VestaObject(double r0_in, double r1_in, double r2_in)
	{
		atom = NULL;
		r[0] = r0_in;
		r[1] = r1_in;
		r[2] = r2_in;
	}

	void Add_Atom(Atom Atom_in)
	{
		Atom* Temp = atom;
		if (atom == NULL)
			atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.e, Atom_in.s);

		else
		{
			while (atom->Next != NULL)
			{
				atom = atom->Next;
			}
			atom_count++;
			atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.e, Atom_in.s);
			atom = Temp;
		}
	}

	void Add_Atom(Atom* Atom_in)
	{
		Atom* Temp = atom;
		if (atom == NULL)
			atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->e, Atom_in->s);

		else
		{
			while (atom->Next != NULL)
			{
				atom = atom->Next;
			}
			atom_count++;
			atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->e, Atom_in->s);
			atom = Temp;
		}
	}

	double PerpLJForce(Atom* Other_Atom)
	{
		Atom* temp = atom;
		double LJForce = 0;

		while (temp != NULL)
		{
			LJForce += temp->PerpLJForce(Other_Atom);
			temp = temp->Next;
		}
		return LJForce;
	}

	double LJForce(Atom* Other_Atom)
	{
		if (Dist(Other_Atom) > 5 * atom->s)
			return(0);

		Atom* temp = atom;
		double LJForce = 0;

		while (temp != NULL)
		{
			LJForce += temp->LJForce(Other_Atom);
			temp = temp->Next;
		}
		return LJForce;
	}

	double LJPot(Atom* Other_Atom)
	{
		Atom* temp = atom;
		double LJPot = 0;

		while (temp != NULL)
		{
			LJPot += temp->LJPot(Other_Atom);
			temp = temp->Next;
		}
		return LJPot;
	}

	double Dist(Atom* Other_Atom)
	{
		return (double)sqrt(pow((Other_Atom->r[2] - r[2]), 2) + pow((Other_Atom->r[1] - r[1]), 2) + pow((Other_Atom->r[0] - r[0]), 2));
	}

	void Print_Atoms()
	{
		Atom* Temp = atom;

		std::ofstream out("VestaObj.dat");
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat

		while (atom != NULL)
		{
			atom->Print();
			atom = atom->Next;
		}

		atom = Temp;
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	void Print()
	{
		cout << r[0] << " " << r[1] << " " << r[2] << "\n";
	}

	void Print(int var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << var << "\n";
	}

	void Print(double var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << var << "\n";
	}

	void Print(double var1, double var2)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << var1 << " " << var2 << "\n";
	}

	void Move(double x_in, double y_in, double z_in)
	{
		Atom* temp = atom;
		while (temp != NULL)
		{
			temp->r[0] += x_in;
			temp->r[1] += y_in;
			temp->r[2] += z_in;
			temp = temp->Next;
		}
		r[0] += x_in;
		r[1] += y_in;
		r[2] += z_in;
	}

	void Import(string filename_in)
	{
		string line;
		string vesta_filename = filename_in;
		string xyz_filename = filename_in;
		string element;
		ifstream fin;
		double rotation[4][4];
		double rotation2[4][4] =
		{
			{ 1,		0,	  0,	 0 },
			{ 0,		0,	 -1,	 0 },
			{ 0,	    1,	  0,	 0 },
			{ 0,		0,	  0,	 0 }
		};
		ResetOrigin();
		double epsilon = 0;
		double sigma = 0;
		int i = 1;
		Atom* Temp = new Atom(0,0,0);
		double initial[3];

		for (int i = 0;i < 3;i++)
		{
			initial[i] = r[i];
			r[i] = 0;
		}



		vesta_filename.erase(0, vesta_filename.find_last_of("\\", vesta_filename.length()) + 1);
		vesta_filename.erase(vesta_filename.find_last_of(".", vesta_filename.length() + 1), vesta_filename.length() + 1);

		xyz_filename.erase(0, xyz_filename.find_last_of("\\", xyz_filename.length()) + 1);
		xyz_filename.erase(xyz_filename.find_last_of(".", xyz_filename.length() + 1), xyz_filename.length() + 1);

		vesta_filename.append(".vesta");
		xyz_filename.append(".xyz");

		fin.open(vesta_filename);

		if (fin.is_open())
		{
			while (getline(fin, line) && line.compare("SCENE") != 0)
			{
			
			}
		}
		else
		{
			cout << vesta_filename << " not found" << endl;
			cin.ignore();
			exit(1);
		}

		for (int i = 0; i < 4; i++)
		{
		 fin >> rotation[i][0]  >> rotation[i][1]  >> rotation[i][2]  >> rotation[i][3];
		}

		fin.close();

		fin.open(xyz_filename);

		if (fin.is_open())
		{
			getline(fin, line);
			getline(fin, line);

			while (fin>>element>>Temp->r[0]>> Temp->r[1]>> Temp->r[2])
			{

				if (element.compare("C") == 0)
				{
					Temp->e = 1.0456*pow(10, -21);	//epsilon	[Joules]	(LJ Parameter) 
					Temp->s = 4.0;					//sigma		[Angstroms] (LJ Parameter) 
				}
				else if (element.compare("Pt") == 0)
				{
					Temp->e = 3.534*pow(10, -21);		//epsilon	[Joules]	(LJ Parameter) 
					Temp->s = 2.95;					//sigma		[Angstroms] (LJ Parameter) 
				}
				else
				{
					cout << "Lennard jones parameters not stored for: " << element << ". Assuming Carbon"<< endl;
					cout << "Call Import(filename, epsilon, sigma)" << endl;
					element = "C";
					Temp->e = 1.0456*pow(10, -21);	//epsilon	[Joules]	(LJ Parameter) 
					Temp->s = 4.0;
				}

				cout << "Atom " << setw(2) << i << ": " << element << endl;
				i++;

				Add_Atom(Temp);
			}
		}
		else
		{
			cout << xyz_filename << " not found" << endl;
			cin.ignore();
			exit(1);
		}
		fin.close();

		Rotate(rotation);
	
		Rotate(rotation2);

		for (int i = 0;i < 3;i++)
			r[i] = initial[i];
		
		ResetOrigin();
		
	}

	void Rotate(double rotation[4][4])
	{
		double r_rot[3];
		Atom* temp = atom;
		while (temp != NULL)
		{
			r_rot[0] = rotation[0][0] * temp->r[0] + rotation[0][1] * temp->r[1] + rotation[0][2] * temp->r[2];
			r_rot[1] = rotation[1][0] * temp->r[0] + rotation[1][1] * temp->r[1] + rotation[1][2] * temp->r[2];
			r_rot[2] = rotation[2][0] * temp->r[0] + rotation[2][1] * temp->r[1] + rotation[2][2] * temp->r[2];

			for (int i = 0;i < 4;i++)
			{
				for (int j = 0;j < 4;j++)
				{
					orientation[i][j] = rotation[i][j]*orientation[j][i];
				}
			}

			temp->r[0] = r_rot[0];
			temp->r[1] = r_rot[1];
			temp->r[2] = r_rot[2];

			temp = temp->Next;
		}
	}

	void Import(string filename_in, double epsilon, double sigma)
	{

		string line;
		string vesta_filename = filename_in;
		string xyz_filename = filename_in;
		string element;
		ifstream fin;
		double orientation[4][4];
		double r_in[3];
		double r_rot[3];


		vesta_filename.erase(0, vesta_filename.find_last_of("\\", vesta_filename.length()) + 1);
		vesta_filename.erase(vesta_filename.find_last_of(".", vesta_filename.length() + 1), vesta_filename.length() + 1);

		xyz_filename.erase(0, xyz_filename.find_last_of("\\", xyz_filename.length()) + 1);
		xyz_filename.erase(xyz_filename.find_last_of(".", xyz_filename.length() + 1), xyz_filename.length() + 1);

		vesta_filename.append(".vesta");
		xyz_filename.append(".xyz");

		fin.open(vesta_filename);

		if (fin.is_open())
		{
			while (getline(fin, line) && line.compare("SCENE") != 0)
			{

			}
		}
		else
		{
			cout << vesta_filename << " not found" << endl;
			cin.ignore();
			exit(1);
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				fin >> orientation[i][j];
			}
		}

		fin.close();

		fin.open(xyz_filename);

		if (fin.is_open())
		{
			getline(fin, line);
			getline(fin, line);

			while (fin >> element >> r_in[0] >> r_in[1] >> r_in[2])
			{
				r_rot[0] = orientation[0][0] * r_in[0] + orientation[0][1] * r_in[1] + orientation[0][2] * r_in[2];
				r_rot[1] = orientation[1][0] * r_in[0] + orientation[1][1] * r_in[1] + orientation[1][2] * r_in[2];
				r_rot[2] = orientation[2][0] * r_in[0] + orientation[2][1] * r_in[1] + orientation[2][2] * r_in[2];


				Atom temp(r_rot[0], r_rot[2], r_rot[1], epsilon, sigma);
				Add_Atom(temp);
			}
		}
		else
		{
			cout << xyz_filename << " not found" << endl;
			cin.ignore();
			exit(1);
		}
		fin.close();

		ResetOrigin();
	}

	void ResetOrigin()
	{
		Atom* temp = atom;
		double r_min[3] = {0, 0, 1000000};
		int i = 0;

		while (atom != NULL)
		{
			r_min[0] += atom->r[0];
			r_min[1] += atom->r[1];
			i++;
			if (atom->r[2] < r_min[2])
			{
				r_min[2] = atom->r[2];
			}

			atom = atom->Next;
		}
		atom = temp;

		r_min[0] = r_min[0] / i;
		r_min[1] = r_min[1] / i;
				
		while (atom != NULL)
		{
			atom->r[0] = atom->r[0] - r_min[0];
			atom->r[1] = atom->r[1] - r_min[1];
			atom->r[2] = atom->r[2] - r_min[2];
			atom = atom->Next;
		}

		r[0] = r_min[0];
		r[1] = r_min[1];
		r[2] = r_min[2];

		atom = temp;
	}

};

class Unit_Cell
{
private:


	
public:

	double r[3];	//position vector of unit cell 
	Unit_Cell* Next;
	Atom* first_atom;

	Unit_Cell()		//Constructor
	{
		first_atom = NULL;
		Next = NULL;
		r[0] = 0.0;
		r[1] = 0.0;
		r[2] = 0.0;
	}

	Unit_Cell(double r0_in, double r1_in, double r2_in)
	{
		first_atom = NULL;
		Next = NULL;
		r[0] = r0_in;
		r[1] = r1_in;
		r[2] = r2_in;
	}
	
	void Add_Atom(Atom Atom_in)
	{
		Atom* Temp = first_atom;
		if (first_atom == NULL)
			first_atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.e, Atom_in.s);

		else
		{
			while (first_atom->Next != NULL)
			{
				first_atom = first_atom->Next;
			}
			first_atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.e, Atom_in.s);
			first_atom = Temp;
		}
	}

	void Add_Atom(Atom* Atom_in)
	{
		Atom* Temp = first_atom;
		if (first_atom == NULL)
			first_atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->e, Atom_in->s);

		else
		{
			while (first_atom->Next != NULL)
			{
				first_atom = first_atom->Next;
			}
			first_atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->e, Atom_in->s);
			first_atom = Temp;
		}
	}

	void Print()
	{
		int atom_number = 0;
		Atom* temp = first_atom;
		while (temp != NULL)
		{
			temp->Print(atom_number);
			temp = temp->Next;
			atom_number++;
		}
	}

	double PerpLJForce(VestaObject* Obj_in)
	{
		Atom* temp = first_atom;
		double Force = 0;

		while (temp != NULL)
		{
			Force += Obj_in->PerpLJForce(temp);
			temp = temp->Next;
		}
		return Force;
	}

	double LJPot(VestaObject* Obj_in)
	{
		Atom* temp = first_atom;
		double pot = 0;

		while (temp != NULL)
		{
			pot += Obj_in->LJPot(temp);
			temp = temp->Next;
		}
		return pot;
	}
	
	double LJForce(VestaObject* Obj_in)
	{	
		
		Atom* temp = first_atom;
		double Force = 0;

		while (temp != NULL)
		{
			Force += Obj_in->LJForce(temp);
			temp = temp->Next;
		}
		return Force;
	}

};

class Surface
{
private:

	Unit_Cell* first_cell;
	double a1[3];	//////////////////////////////////////////////////////////
	double a2[3];	//		Lattice vectors a1, a2 and a3 follow from		//
	double a3[3];	//	crystallographic convention and have units of nm	//
	int a1_cells;	//////////////////////////////////////////////////////////
	int a2_cells;
	int a3_cells;

	void Tile_Space()
	{
		Unit_Cell* temp_cell;
		for (int k = 0; k < a3_cells; k++)
		{
			for (int j = 0; j < a2_cells; j++)
			{
				for (int i = 0; i < a1_cells; i++)
				{
					if (i == 0 && j == 0 && k == 0)
					{
						i++;
					}
					temp_cell = first_cell;
					first_cell = new Unit_Cell((a1[0] * i + a2[0] * j + a3[0] * (k % 2)), (a1[1] * i + a2[1] * j + a3[1] * (k % 2)), (a1[2] * i + a2[2] * j + a3[2] * k));
					cell_count++;
					first_cell->Next = temp_cell;
				}
			}
		}
	}

	void Fill_Tiled_Space(Unit_Cell* Cell_in)
	{
		Unit_Cell* temp_cell = first_cell;
		Atom* temp_atom = Cell_in->first_atom;
		double a1_factor;
		double a2_factor;
		double a3_factor;

		while (Cell_in->first_atom != NULL)
		{
			a1_factor = Cell_in->first_atom->r[0];
			a2_factor = Cell_in->first_atom->r[1];
			a3_factor = Cell_in->first_atom->r[2];
			Cell_in->first_atom->r[0] = (a1_factor * a1[0]) + (a2_factor * a2[0]) + (a3_factor * a3[0]);
			Cell_in->first_atom->r[1] = (a1_factor * a1[1]) + (a2_factor * a2[1]) + (a3_factor * a3[1]);
			Cell_in->first_atom->r[2] = (a1_factor * a1[2]) + (a2_factor * a2[2]) + (a3_factor * a3[2]);
			Cell_in->first_atom = Cell_in->first_atom->Next;
		}


		Cell_in->first_atom = temp_atom;

		while (first_cell->Next != NULL)
		{
			while (Cell_in->first_atom != NULL)
			{
				first_cell->Add_Atom(Cell_in->first_atom);
				atom_count++;
				Cell_in->first_atom = Cell_in->first_atom->Next;
			}
			Cell_in->first_atom = temp_atom;
			first_cell = first_cell->Next;
		}
		first_cell = temp_cell;
	}

	void TipHeightMap(VestaObject* Obj_in, double XMax, double XMin, double YMax, double YMin, double ZIni, double Setpoint, double ZStep, double FRes)
	{
		int progress = 0;
		double Prev = 0;
		double Cur = 0;
		double xstep = (XMax - XMin) / 250;
		double ystep = (YMax - YMin) / 250;

		if (XMax == XMin)
			xstep = 1;

		if (YMax == YMin)
			ystep = 1;

		double original_position[3] = { Obj_in->r[0], Obj_in->r[1], Obj_in->r[2] };

		cout << "x [Å] " << "y [Å] " << "z [Å] " << "dz " << endl;

		Obj_in->Move(XMin - Obj_in->r[0], YMin - Obj_in->r[1], ZIni - Obj_in->r[2]);


		while (Obj_in->r[0] <= XMax)
		{

			Prev = TipHeightCalc(Obj_in, Setpoint, ZStep, FRes);
			Cur = TipHeightCalc(Obj_in, Setpoint, ZStep, FRes);
			while (Obj_in->r[1] <= YMax)
			{
				Cur  = TipHeightCalc(Obj_in, Setpoint, ZStep, FRes);
				Obj_in->Print(100*(Cur-Prev));
				Prev = Cur;
				Obj_in->Move(0, ystep, 0);
			}
			Obj_in->Move(xstep, YMin - Obj_in->r[1], 0);

			progress = (int)((Obj_in->r[0] - XMin) * 100 / (XMax - XMin));
			printf("Topology %d%% Complete\r", progress);
			if (YMax != YMin)
			{
				cout << "\n";
			}
		}
		printf("\n");
		Obj_in->Move(original_position[0] - Obj_in->r[0], original_position[1] - Obj_in->r[1], original_position[2] - Obj_in->r[2]);
	}

	double TipHeightCalc(VestaObject* Obj_in, double Setpoint, double ZStep, double FRes)
	{
		double Error = LJForce(Obj_in)-Setpoint;
		bool check1 = 0;
		bool check2 = 0;

		while (abs(Error) >  FRes)
		{
			while ((check1 * check2) != 1)
			{
				if (Error > 0)
				{
					Obj_in->Move(0, 0, ZStep);
					check1 = 1;
				}

				else
				{
					Obj_in->Move(0, 0, -1.0*ZStep);
					check2 = 1;
				}
				Error = LJForce(Obj_in) - Setpoint;
			}

		    	ZStep = ZStep /5;
				check1 = 0;
				check2 = 0;
		}
		return(Obj_in->r[2]);
	}

	double LJPot(VestaObject* Obj_in)
	{
		Unit_Cell* temp = first_cell;
		double LJPot = 0;

		while (temp != NULL)
		{
			LJPot += temp->LJPot(Obj_in);
			temp = temp->Next;
		}
		return LJPot;
	}

	double LJForce(VestaObject* Obj_in)
	{
	

		Unit_Cell* temp = first_cell;
		double Force = 0;

		while (temp != NULL)
		{
			Force += temp->LJForce(Obj_in);
			temp = temp->Next;
		}
		return Force;
	}

	double PerpLJForce(VestaObject* Obj_in)
	{
		Unit_Cell* temp = first_cell;
		double Force = 0;

		while (temp != NULL)
		{
			Force += temp->PerpLJForce(Obj_in);
			temp = temp->Next;
		}
		return Force;
	}

public:

	int cell_count = 0;
	int atom_count = 0;

	Surface(double a1_in[3], double a2_in[3], double a3_in[3], int a1_cells_in, int a2_cells_in, int a3_cells_in, Unit_Cell* Cell_in)
	{
		a1[0] = a1_in[0];
		a1[1] = a1_in[1];
		a1[2] = a1_in[2];

		a2[0] = a2_in[0];
		a2[1] = a2_in[1];
		a2[2] = a2_in[2];

		a3[0] = a3_in[0];
		a3[1] = a3_in[1];
		a3[2] = a3_in[2];

		a1_cells = a1_cells_in;
		a2_cells = a2_cells_in;
		a3_cells = a3_cells_in;

		first_cell = Cell_in;

		Tile_Space();
		Fill_Tiled_Space(Cell_in);

	}

	void Print()
	{
		std::ofstream out("Locations.dat");
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat

		Unit_Cell* temp = first_cell;
		while (temp != NULL)
		{
			temp->Print();
			temp = temp->Next;
		}
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}
	 
	void TipHeight(VestaObject* Obj_in, double Setpoint, double ZRes, double FRes, double Area)
	{
		double XMax = round(((a1[0] * a1_cells + a2[0] * a2_cells) / 2 + Area * (abs(a1[0]) + abs(a2[0]))));
		double XMin = round(((a1[0] * a1_cells + a2[0] * a2_cells) / 2 - Area * (abs(a1[0]) + abs(a2[0]))));
		double YMax = round(((a1[1] * a1_cells + a2[1] * a2_cells) / 2 + Area * (abs(a1[1]) + abs(a2[1]))));
		double YMin = round(((a1[1] * a1_cells + a2[1] * a2_cells) / 2 - Area * (abs(a1[1]) + abs(a2[1]))));
		double Z = Obj_in->r[2];
		std::stringstream buffer;
		std::ofstream out("Topology.dat");

		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(buffer.rdbuf());

		TipHeightMap(Obj_in, XMax, XMin, YMax, YMin, Z, Setpoint, ZRes, FRes);

		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat
		cout << buffer.str() << endl;
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	void SurfaceForce(VestaObject* Obj_in, double Area)
	{
		double XMax = round((a1[0] * a1_cells + a2[0] * a2_cells) / 2 + Area * (abs(a1[0]) + abs(a2[0])));
		double XMin = round((a1[0] * a1_cells + a2[0] * a2_cells) / 2 - Area * (abs(a1[0]) + abs(a2[0])));
		double YMax = round((a1[1] * a1_cells + a2[1] * a2_cells) / 2 + Area * (abs(a1[1]) + abs(a2[1])));
		double YMin = round((a1[1] * a1_cells + a2[1] * a2_cells) / 2 - Area * (abs(a1[1]) + abs(a2[1])));
		double Z = Obj_in->r[2];
		std::stringstream buffer;
		std::ofstream out("Surface.dat");

		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(buffer.rdbuf());

		XYZForce(Obj_in, XMax, XMin, YMax, YMin, Z, Z, 1);
				
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat
		cout << buffer.str() << endl;
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	void ForceCurve(VestaObject* Obj_in, double ZMin, double ZMax)
	{

		double XMax = round((((a1_cells / 2) + 0.75)*a1[0] + ((a2_cells / 2) + 0.75)*a2[0]));
		double XMin = round((((a1_cells / 2) - 1.25)*a1[0] + ((a2_cells / 2) - 1.25)*a2[0]));
		double Y = 0;
		std::stringstream buffer;
		std::ofstream out("Force_Curve.dat");

		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(buffer.rdbuf());

		XYZForce(Obj_in, XMax, XMin, Y, Y, ZMax, ZMin, 1);

		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat
		cout << buffer.str() << endl;
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}
	
	void XYZForce(VestaObject* Obj_in, double XMax, double XMin, double YMax, double YMin, double ZMax, double ZMin, bool Corrugation)
	{
		double Average = 0;
		int progress = 0;
		int i = 0;
		double xstep = (XMax - XMin) / 250;
		double ystep = (YMax - YMin) / 250;
		double zstep = (ZMax - ZMin) / 250;

		if (XMax == XMin)
			xstep = 1;
		
		if (YMax == YMin)
			ystep = 1;

		if (ZMax == ZMin)
			zstep = 1;
		
		double original_position[3] = { Obj_in->r[0], Obj_in->r[1], Obj_in->r[2] };

		cout << "x [Å] " << "y [Å] " << "z [Å] " << "U [nN] " << endl;

		Obj_in->Move(XMin - Obj_in->r[0], YMin - Obj_in->r[1], ZMin - Obj_in->r[2]);

		while (Obj_in->r[2] <= ZMax)
		{
			while (Obj_in->r[0] <= XMax && Corrugation==1)
			{
				while (Obj_in->r[1] <= YMax)
				{
					Average += LJForce(Obj_in);
					i++;
					Obj_in->Move(0, ystep*25, 0);
				}
				Obj_in->Move(xstep*25, YMin- Obj_in->r[1], 0);
			}
			
			Average = Average / i;
			Obj_in->Move(XMin- Obj_in->r[0], YMin - Obj_in->r[1], 0);

			while (Obj_in->r[0] <= XMax)
			{
				while (Obj_in->r[1] <= YMax)
				{
					Obj_in->Print(LJForce(Obj_in)-Average);
					Obj_in->Move(0, ystep, 0);
				}
				Obj_in->Move(xstep, YMin - Obj_in->r[1], 0);
				if (YMax != YMin)
				{
					cout << "\n";
				}

				if (ZMax < ZMin*1.01)
				{
					progress = (int)((Obj_in->r[0] - XMin) * 100 / (XMax - XMin));
					printf("Force calculation %d%% Complete\r", progress);
				}
			}

			if (ZMax > ZMin*1.01)
			{
				progress = (int)((Obj_in->r[2] - ZMin) * 100 / (ZMax - ZMin));
				printf("Force calculation %d%% Complete\r", progress);
			}

			Obj_in->Move(XMin- Obj_in->r[0], YMin - Obj_in->r[1], zstep);
			cout << "\n";
			Average = 0;
			i = 0;
		}
		printf("\n");
		Obj_in->Move(original_position[0] - Obj_in->r[0], original_position[1] - Obj_in->r[1], original_position[2] - Obj_in->r[2]);
	}
	
	void Add_Defect(Unit_Cell* Defect_Cell)
	{
		Unit_Cell* Temp=first_cell;
		while (first_cell->Next != NULL)
			first_cell = first_cell->Next;

		first_cell->Next = Defect_Cell;

		first_cell = Temp;
	}

	void Add_Defect(VestaObject* Object_In)
	{
		Atom* TempAtom = Object_In->atom;
		Unit_Cell* DefectCell = new Unit_Cell(0, 0, 0);


		while (TempAtom != NULL)
		{
			DefectCell->Add_Atom(TempAtom);
			TempAtom = TempAtom->Next;
		}
		Add_Defect(DefectCell);
	}
};

#endif
