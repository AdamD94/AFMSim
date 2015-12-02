#ifndef Class_Definitions
#define Class_Definitions

class Atom
{
public:
	double r[3];
	double	e;
	double	s;
	double temp;
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
		if (temp > 3 * s)
			return 0;

		return (double)	 4 * e*(pow((s /temp), 12) - pow((s / temp), 6));
	}

	double LJForce(Atom* Other_Atom)
	{
		temp = Dist(Other_Atom);
		if (temp > 3 * s)
			return 0;

		return (double) (4.0 / temp)*e*(12 * pow((s / temp), 12) - 6 * pow((s / temp), 6))*pow(10, 19);
	}
	
	double PerpLJForce(Atom* Other_Atom)
	{
		return (double) LJForce(Other_Atom)*(abs(Other_Atom->r[2]-r[2]))/(Dist(Other_Atom));
	}

};

class Tip
{
public:
	Atom* atom;
	double r[3];

	Tip(double r0_in, double r1_in, double r2_in)
	{
		atom = NULL;
		r[0] = r0_in;
		r[1] = r1_in;
		r[2] = r2_in;
	}

	double PerpLJForce(Atom* Other_Atom)
	{
		Atom* Temp = atom;
		double LJForce = 0;

		while (atom != NULL)
		{
			LJForce += atom->PerpLJForce(Other_Atom);
			atom = atom->Next;
		}
		atom = Temp;
		return LJForce;
	}

	double LJForce(Atom* Other_Atom)
	{
		Atom* Temp = atom;
		double LJForce = 0;

		while (atom != NULL)
		{
			LJForce += atom->LJForce(Other_Atom);
			atom = atom->Next;
		}
		atom = Temp;
		return LJForce;
	}

	double LJPot(Atom* Other_Atom)
	{
		Atom* Temp = atom;
		double LJPot = 0;

		while (atom != NULL)
		{
			LJPot += atom->LJPot(Other_Atom);
			atom = atom->Next;
		}
		atom = Temp;
		return LJPot;
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
			atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->e, Atom_in->s);
			atom = Temp;
		}
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

	void MoveTip(double x_in, double y_in, double z_in)
	{
		Atom* Temp = atom;
		while (atom != NULL)
		{
			atom->r[0] += x_in;
			atom->r[1] += y_in;
			atom->r[2] += z_in;
			atom = atom->Next;
		}
		atom = Temp;
		r[0] += x_in;
		r[1] += y_in;
		r[2] += z_in;
	}

};

class Unit_Cell
{
public:
	Atom* first_atom;
	Unit_Cell* Next;
	double r[3];	//position vector of unit cell 

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
		Atom* Temp = first_atom;
		while (first_atom != NULL)
		{
			first_atom->Print(atom_number);
			first_atom = first_atom->Next;
			atom_number++;
		}
		first_atom = Temp;
	}

	double PerpLJForce(Tip* Tip_in)
	{
		Atom* Temp = first_atom;
		double Force = 0;

		while (first_atom != NULL)
		{
			Force += Tip_in->PerpLJForce(first_atom);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return Force;
	}

	double LJPot(Tip* Tip_in)
	{

		Atom* Temp = first_atom;
		double pot = 0;

		while (first_atom != NULL)
		{
			pot += Tip_in->LJPot(first_atom);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return pot;
	}
	
	double LJForce(Tip* Tip_in)
	{
		Atom* Temp = first_atom;
		double Force = 0;

		while (first_atom != NULL)
		{
			Force += Tip_in->LJForce(first_atom);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return Force;
	}

};

class Surface
{
public:

	Unit_Cell* first_cell;

	double a1[3];	//////////////////////////////////////////////////////////
	double a2[3];	//		Lattice vectors a1, a2 and a3 follow from		//
	double a3[3];	//	crystallographic convention and have units of nm	//
	int a1_cells;	//////////////////////////////////////////////////////////
	int a2_cells;
	int a3_cells;

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
				Cell_in->first_atom = Cell_in->first_atom->Next;
			}
			Cell_in->first_atom = temp_atom;
			first_cell = first_cell->Next;
		}
		first_cell = temp_cell;
	}

	void Print()
	{
		std::ofstream out("Locations.dat");
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat

		Unit_Cell* temp_cell = first_cell;
		while (first_cell != NULL)
		{
			first_cell->Print();
			first_cell = first_cell->Next;
		}
		first_cell = temp_cell;
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	double LJPot(Tip* Tip_in)
	{
		Unit_Cell* Temp = first_cell;
		double LJPot = 0;

		while (first_cell != NULL)
		{
			LJPot += first_cell->LJPot(Tip_in);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
		return LJPot;
	}

	double LJForce(Tip* Tip_in)
	{
		Unit_Cell* Temp = first_cell;
		double Force = 0;

		while (first_cell != NULL)
		{
			Force += first_cell->LJForce(Tip_in);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
		return Force;
	}

	double PerpLJForce(Tip* Tip_in)
	{
		Unit_Cell* Temp = first_cell;
		double Force = 0;

		while (first_cell != NULL)
		{
			Force += first_cell->PerpLJForce(Tip_in);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
		return Force;
	}

	void SurfaceForce(Tip* Tip_in)
	{
		double Average = 0;
		double XMax = (a1[0] * a1_cells + a2[0] * a2_cells) / 2 + 1 * (abs(a1[0]) + abs(a2[0]));
		double XMin = (a1[0] * a1_cells + a2[0] * a2_cells) / 2 - 1 * (abs(a1[0]) + abs(a2[0]));
		double YMax = (a1[1] * a1_cells + a2[1] * a2_cells) / 2 + 1 * (abs(a1[1]) + abs(a2[1]));
		double YMin = (a1[1] * a1_cells + a2[1] * a2_cells) / 2 - 1 * (abs(a1[1]) + abs(a2[1]));
		double Z = Tip_in->r[2];

		std::ofstream out("Surface.dat");
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat

		XYZForce(Tip_in, XMax, XMin, YMax, YMin, Z, Z);

		cout << endl;
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	void ForceCurve(Tip* Tip_in)
	{

		double Average = 0;
		double XMax = (((a1_cells / 2) + 0.75)*a1[0] + ((a2_cells / 2) + 0.75)*a2[0]);
		double XMin = (((a1_cells / 2) - 1.25)*a1[0] + ((a2_cells / 2) - 1.25)*a2[0]);
		double Y = 0;
		double ZMax = 4.0;
		double ZMin = 3.0;

		std::ofstream out("Force_Curve.dat");
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat

		XYZForce(Tip_in, XMax, XMin, Y, Y, ZMax, ZMin);

		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	void XYZForce(Tip* Tip_in, double XMax, double XMin, double YMax, double YMin, double ZMax, double ZMin)
	{
		double Average = 0;

		int i = 0;
		double xstep = (XMax - XMin) / 750;
		double ystep = (YMax - YMin) / 750;
		double zstep = (ZMax - ZMin) / 750;

		if (XMax == XMin)
			xstep = 1;
		
		if (YMax == YMin)
			ystep = 1;

		if (ZMax == ZMin)
			zstep = 1;
		
		double original_position[3] = { Tip_in->r[0], Tip_in->r[1], Tip_in->r[2] };

		cout << "x [Å] " << "y [Å] " << "z [Å] " << "U [nN] " << endl;

		Tip_in->MoveTip(XMin - Tip_in->r[0], YMin - Tip_in->r[1], ZMin - Tip_in->r[2]);

		while (Tip_in->r[2] <= ZMax)
		{
			while (Tip_in->r[0] <= XMax)
			{
				while (Tip_in->r[1] <= YMax)
				{
					Average += LJForce(Tip_in);
					i++;
					Tip_in->MoveTip(0, ystep*25, 0);
				}
				Tip_in->MoveTip(xstep*25, YMin-Tip_in->r[1], 0);
			}
			
			Average = Average / i;
			Tip_in->MoveTip(XMin-Tip_in->r[0], YMin - Tip_in->r[1], 0);

			while (Tip_in->r[0] <= XMax)
			{
				while (Tip_in->r[1] <= YMax)
				{
					Tip_in->Print(LJForce(Tip_in)-Average);
					Tip_in->MoveTip(0, ystep, 0);
				}
				Tip_in->MoveTip(xstep, YMin - Tip_in->r[1], 0);
				if (YMax != YMin)
				{
					cout << "\n";
				}
			}
			Tip_in->MoveTip(XMin-Tip_in->r[0], YMin - Tip_in->r[1], zstep);
			cout << "\n";
			Average = 0;
			i = 0;
		}

		
		Tip_in->MoveTip(original_position[0] - Tip_in->r[0], original_position[1] - Tip_in->r[1], original_position[2] - Tip_in->r[2]);
	}

	void PerpSurfaceForce(Tip* Tip_in)
	{
		double Average = 0;
		double XMax = (a1[0] * a1_cells + a2[0] * a2_cells) / 2 + 1 * (abs(a1[0]) + abs(a2[0]));
		double XMin = (a1[0] * a1_cells + a2[0] * a2_cells) / 2 - 1 * (abs(a1[0]) + abs(a2[0]));
		double YMax = (a1[1] * a1_cells + a2[1] * a2_cells) / 2 + 1 * (abs(a1[1]) + abs(a2[1]));
		double YMin = (a1[1] * a1_cells + a2[1] * a2_cells) / 2 - 1 * (abs(a1[1]) + abs(a2[1]));
		int k = 0;
		double step = (XMax - XMin) / 750;
		double original_position[3] = { Tip_in->r[0], Tip_in->r[1], Tip_in->r[2] };

		std::ofstream out("Surface.dat");
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to AFM.dat

		cout << "x " << "y " << "z " << "U " << endl;

		Tip_in->MoveTip(XMin - Tip_in->r[0], YMin - Tip_in->r[1], 0);

		Tip_in->MoveTip(0, YMin - Tip_in->r[1], 0);
		while (Tip_in->r[1] < YMax)
		{
			Tip_in->MoveTip(XMin - Tip_in->r[0], 0, 0);
			while (Tip_in->r[0] < XMax)
			{
				Average += PerpLJForce(Tip_in);
				Tip_in->MoveTip(step * 10, 0, 0);
				k++;
			}
			Tip_in->MoveTip(0, step * 10, 0);
		}

		Average = Average / k;
		Tip_in->MoveTip(0, YMin - Tip_in->r[1], 0);

		while (Tip_in->r[1] < YMax)
		{
			Tip_in->MoveTip(XMin - Tip_in->r[0], 0, 0);
			while (Tip_in->r[0] < XMax)
			{
				Tip_in->Print(PerpLJForce(Tip_in));
				Tip_in->MoveTip(step, 0, 0);
			}
			Tip_in->MoveTip(0, step, 0);
		}

		cout << endl;
		std::cout.rdbuf(coutbuf); //reset to standard output again
		Tip_in->MoveTip(original_position[0] - Tip_in->r[0], original_position[1] - Tip_in->r[1], original_position[2] - Tip_in->r[2]);
	}

	void PerpForceCurve(Tip* Tip_in)
	{
		double Average = 0;
		double XMax = (((a1_cells / 2) + 0.75)*a1[0] + ((a2_cells / 2) + 0.75)*a2[0]);
		double XMin = (((a1_cells / 2) - 1.25)*a1[0] + ((a2_cells / 2) - 1.25)*a2[0]);
		double YMax = abs(a1[1]) + abs(a2[1]);
		double YMin = 0;
		double ZMax = 4;
		double ZMin = 3;
		int k = 0;
		double zstep = (ZMax - ZMin) / 750;
		double xystep = (XMax - XMin) / 750;
		double original_position[3] = { Tip_in->r[0], Tip_in->r[1], Tip_in->r[2] };

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
					Average += PerpLJForce(Tip_in);
					k++;
					Tip_in->MoveTip(xystep * 10, 0, 0);
				}
				Tip_in->MoveTip(0, xystep * 10, 0);
			}
			Average = Average / k;

			Tip_in->MoveTip(XMin - Tip_in->r[0], YMin - Tip_in->r[1], 0);

			while (Tip_in->r[0] < XMax)
			{

				Tip_in->Print(PerpLJForce(Tip_in) - Average);
				Tip_in->MoveTip(xystep, 0, 0);
			}

			Average = 0;
			k = 0;
			Tip_in->MoveTip(XMin - Tip_in->r[0], 0, zstep);
		}

		cout << endl;
		std::cout.rdbuf(coutbuf); //reset to standard output again
		Tip_in->MoveTip(original_position[0] - Tip_in->r[0], original_position[1] - Tip_in->r[1], original_position[2] - Tip_in->r[2]);
	}

};

#endif