#ifndef Class_Definitions
#define Class_Definitions

class Atom
{
public:
	double r[3];
	double radius;
	double m;
	double epsilon;
	double sigma;
	Atom* Next;

	Atom(double x_in, double y_in, double z_in, double radius_in, double m_in, double epsilon_in, double sigma_in)
	{
		r[0] = x_in;
		r[1] = y_in;
		r[2] = z_in;
		radius = radius_in;
		m = m_in;
		epsilon = epsilon_in;
		sigma = sigma_in;
		Next = NULL;
	}

	void Print()
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius <<" "<< m << endl;
	}

	void Print(double force)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius << " " << m << " " << force << endl;
	}

	void Print(int atom)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius << " " << m <<" "<<atom<< endl;
	}

	double LJForce(Atom* Other_Atom)
	{
		return(0);
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
		if (first_atom == NULL)
			first_atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.radius, Atom_in.m, Atom_in.epsilon, Atom_in.sigma);

		else
			first_atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.radius, Atom_in.m, Atom_in.epsilon, Atom_in.sigma);
	}

	void Add_Atom(Atom* Atom_in)
	{
		if (first_atom == NULL)
			first_atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m, Atom_in->epsilon, Atom_in->sigma);

		else
			first_atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m, Atom_in->epsilon, Atom_in->sigma);
	}

	void Print()
	{
		int atom_number = 0;
		Atom* temp_atom = first_atom;
		while (first_atom != NULL)
		{
			first_atom->Print(atom_number);
			first_atom = first_atom->Next;
			atom_number++;
		}
		first_atom = temp_atom;
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
		for (double k = 0; k < a3_cells; k++)
		{
			for (double j = 0; j < a2_cells; j++)
			{
				for (double i = 0; i < a1_cells; i++)
				{
					if (i == 0 && j == 0 && k == 0)
					{ 
						i++;
					}
					temp_cell = first_cell;
					first_cell= new Unit_Cell((a1[0]*i + a2[0]*j + a3[0] * k), (a1[1] * i + a2[1] * j + a3[1] * k) , (a1[2] * i + a2[2] * j + a3[2] * k));
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



};

class Surface
{
public:
	Lattice* first_Lattice;
	Atom* first_Atom;
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
};

#endif