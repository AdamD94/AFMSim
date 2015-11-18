#ifndef Class_Definitions
#define Class_Definitions

class Atom
{
public:
	double r[3];
	double radius;
	double m;
	char elem[6];
	Atom* Next;

	Atom(double x_in, double y_in, double z_in, double radius_in, double m_in, char elem_in[6] )
	{
		r[0] = x_in;
		r[1] = y_in;
		r[2] = z_in;
		radius = radius_in;
		m = m_in;
		strcpy_s(elem,elem_in);
		Next = NULL;
	}

	void Print()
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius <<" "<< m << "\n";
	}

	void Print(double var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius << " " << m << " " << var << "\n";
	}

	void Print(int var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius << " " << m << " " << var << "\n";
	}

	void XYPrint(double var)
	{
		cout << r[0] << " " << r[1] << " " << var << "\n";
	}

	void XYPrint(double var1, double var2)
	{
		cout << r[0] << " " << r[1] << " " << var1 << " " << var2 << "\n";
	}

	double Dist(Atom* Other_Atom)
	{
		double dist = sqrt(pow((Other_Atom->r[2] - r[2]), 2) + pow((Other_Atom->r[1] - r[1]), 2) + pow((Other_Atom->r[0] - r[0]), 2));
		return dist;
	}

	double LJPot(Atom* Other_Atom)
	{
		double e = 3.22*pow(10, -20);
		double s = 0.241;

		double r = Dist(Other_Atom);

		double LJPotential = 4 * e*(pow((s / r), 12) - pow((s / r), 6));

		return LJPotential;
	}

	double LJForce(Atom* Other_Atom)
	{
		double e = 3.534*pow(10, -21);
		double s = 0.295;

		double r = Dist(Other_Atom);

		double LJForce = (-4/r)*e*(12*pow((s / r), 12) - 6*pow((s / r), 6))*pow(10, 18);

		return LJForce;
	}
	
	double PerpLJForce(Atom* Other_Atom)
	{
		double Force = LJForce(Other_Atom)*(Other_Atom->r[2]-r[2])/(Dist(Other_Atom));
		return Force;
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

	void Add_Atom(Atom Atom_in)
	{
		Atom* Temp = atom;
		if (atom == NULL)
			atom = new Atom(Atom_in.r[0], Atom_in.r[1], Atom_in.r[2], Atom_in.radius, Atom_in.m, Atom_in.elem);

		else
		{
			while (atom->Next != NULL)
			{
				atom = atom->Next;
			}

			atom->Next = new Atom(Atom_in.r[0], Atom_in.r[1], Atom_in.r[2], Atom_in.radius, Atom_in.m, Atom_in.elem);
			atom = Temp;
		}
	}

	void Add_Atom(Atom* Atom_in)
	{
		Atom* Temp = atom;
		if (atom == NULL)
			atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m, Atom_in->elem);

		else
		{
			while (atom->Next != NULL)
			{
				atom = atom->Next;
			}
			atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m, Atom_in->elem);
			atom = Temp;
		}
	}

	void Print()
	{
		cout << r[0] << " " << r[1] << " " << r[2] << "\n";
	}

	void Print(int var)
	{
		cout << r[0] << " " << r[1] << " " << r[2] <<" " << var << "\n";
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
	class Lattice;
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
			first_atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.radius, Atom_in.m, Atom_in.elem);

		else
		{
			while (first_atom->Next != NULL)
			{
				first_atom = first_atom->Next;
			}
			first_atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.radius, Atom_in.m, Atom_in.elem);
			first_atom = Temp;
		}
	}

	void Add_Atom(Atom* Atom_in)
	{
		Atom* Temp = first_atom;
		if (first_atom == NULL)
			first_atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m, Atom_in->elem);

		else
		{
			while (first_atom->Next != NULL)
			{
				first_atom = first_atom->Next;
			}
			first_atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m, Atom_in->elem);
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

	void Print(double var)
	{
		Atom* Temp = first_atom;
		while (first_atom->Next != NULL)
		{
			first_atom->Print(var);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
	}

	double LJPot(Atom* atom_in)
	{
		Atom* Temp = first_atom;
		double pot = 0;

		while (first_atom->Next != NULL)
		{
			pot += first_atom->LJPot(atom_in);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return pot;
	}

	double LJForce(Atom* atom_in)
	{
		Atom* Temp = first_atom;
		double Force = 0;

		while (first_atom->Next != NULL)
		{
			Force += first_atom->LJForce(atom_in);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return Force;
	}

	double LJForce(Tip* Tip_in)
	{
		Atom* Temp = first_atom;
		double Force = 0;

		while (first_atom->Next != NULL)
		{
			Force += Tip_in->LJForce(first_atom);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return Force;
	}

	double PerpLJForce(Atom* atom_in)
	{
		Atom* Temp = first_atom;
		double Force = 0;

		while (first_atom->Next != NULL)
		{
			Force += first_atom->PerpLJForce(atom_in);
			first_atom = first_atom->Next;
		}
		first_atom = Temp;
		return Force;
	}

};

class Lattice
{

public:

	Unit_Cell* first_cell;
	Lattice* Next;

	double a1[3];	//////////////////////////////////////////////////////////
	double a2[3];	//		Lattice vectors a1, a2 and a3 follow from		//
	double a3[3];	//	crystallographic convention and have units of nm	//
					//////////////////////////////////////////////////////////

	Lattice(double a1_in[3], double a2_in[3], double a3_in[3])	//Constructor
	{
		first_cell = NULL;
		Next = NULL;

		a1[0] = a1_in[0];
		a1[1] = a1_in[1];
		a1[2] = a1_in[2];

		a2[0] = a2_in[0];
		a2[1] = a2_in[1];
		a2[2] = a2_in[2];

		a3[0] = a3_in[0];
		a3[1] = a3_in[1];
		a3[2] = a3_in[2];
	}
		
	void Tile_Space(int a1_cells, int a2_cells, int a3_cells)
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
					first_cell= new Unit_Cell((a1[0]*i + a2[0]*j + a3[0] * (k%2)), (a1[1] * i + a2[1] * j + a3[1] * (k%2) ) , (a1[2] * i + a2[2] * j + a3[2] * k));
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
		Unit_Cell* temp_cell= first_cell;
		while (first_cell != NULL)
		{
			first_cell->Print();
			first_cell = first_cell->Next;
		}
		first_cell = temp_cell;
	}

	void Print(double var)
	{
		Unit_Cell* Temp = first_cell;

		while (first_cell != NULL)
		{
			first_cell->Print(var);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
	}

	double LJPot(Atom* atom_in)
	{
		Unit_Cell* Temp = first_cell;
		double pot = 0;

		while (first_cell != NULL)
		{
			pot += first_cell->LJPot(atom_in);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
		return pot;
	}

	double LJForce(Atom* atom_in)
	{
		Unit_Cell* Temp = first_cell;
		double Force = 0;

		while (first_cell != NULL)
		{
			Force += first_cell->LJForce(atom_in);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
		return Force;
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

	double PerpLJForce(Atom* atom_in)
	{
		Unit_Cell* Temp = first_cell;
		double Force = 0;

		while (first_cell != NULL)
		{
			Force += first_cell->PerpLJForce(atom_in);
			first_cell = first_cell->Next;
		}
		first_cell = Temp;
		return Force;
	}


};

class Surface
{
public:
	Lattice* first_Lattice;
	Unit_Cell* first_cell;
	Atom* first_Atom;
	Atom* Linked_Atoms;
	int layer;

	Surface()		//Constructor
	{
		first_Lattice = NULL;
		first_Atom = NULL;
	}

	Surface(Lattice* Lattice_in, Unit_Cell* Cell_in, int a1_cells, int	a2_cells, int a3_cells)
	{
		first_Lattice = Lattice_in;
		first_Lattice->first_cell = Cell_in;

		first_Lattice->Tile_Space(a1_cells, a2_cells, a3_cells);
		first_Lattice->Fill_Tiled_Space(Cell_in);
		first_cell = first_Lattice->first_cell;
		first_Atom = first_Lattice->first_cell->first_atom;
	}

	void Print()
	{
		Lattice* temp_Lattice = first_Lattice;
		while (first_Lattice != NULL)
		{
			first_Lattice->Print();
			first_Lattice = first_Lattice->Next;
		}
		first_Lattice = temp_Lattice;
	}

	void Print(double var)
	{
		Lattice* temp_Lattice = first_Lattice;
		while (first_Lattice != NULL)
		{
			first_Lattice->Print(var);
			first_Lattice = first_Lattice->Next;
		}
		first_Lattice = temp_Lattice;
	}

	double LJPot(Atom* atom_in)
	{
		Lattice* temp_Lattice = first_Lattice;
		double pot = 0;
		while (first_Lattice != NULL)
		{
			pot += first_Lattice->LJPot(atom_in);
			first_Lattice = first_Lattice->Next;
		}
		first_Lattice = temp_Lattice;
		return pot;
	}

	double LJForce(Atom* atom_in)
	{
		Lattice* temp_Lattice = first_Lattice;
		double Force = 0;
		while (first_Lattice != NULL)
		{
			Force += first_Lattice->LJForce(atom_in);
			first_Lattice = first_Lattice->Next;
		}
		first_Lattice = temp_Lattice;
		return Force;
	}

	double PerpLJForce(Atom* atom_in)
	{
		Lattice* temp_Lattice = first_Lattice;
		double Force = 0;
		while (first_Lattice != NULL)
		{
			Force += first_Lattice->PerpLJForce(atom_in);
			first_Lattice = first_Lattice->Next;
		}
		first_Lattice = temp_Lattice;
		return Force;
	}

	double LJForce(Tip* Tip_in)
	{
		Lattice* temp_Lattice = first_Lattice;
		double Force = 0;
		while (first_Lattice != NULL)
		{
			Force += first_Lattice->LJForce(Tip_in);
			first_Lattice = first_Lattice->Next;
		}
		first_Lattice = temp_Lattice;
		return Force;
	}

};

#endif