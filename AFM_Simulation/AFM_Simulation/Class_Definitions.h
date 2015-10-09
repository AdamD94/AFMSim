#ifndef Class_Definitions
#define Class_Definitions

class Atom
{
public:
	double r[3];
	double radius;
	double m;
	Atom* Next;

	Atom(double x_in, double y_in, double z_in, double radius_in, double m_in)
	{
		r[0] = x_in;
		r[1] = y_in;
		r[2] = z_in;
		radius = radius_in;
		m = m_in;
		Next = NULL;
	}

	void Print()
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius <<" "<< m << endl;
	}

	void Print(int lattice, int atom)
	{
		cout << r[0] << " " << r[1] << " " << r[2] << " " << radius << " " << m << " " <<lattice<<" "<<atom<< endl;
	}

};

class Unit_Cell
{
public:
	class Lattice;
	Atom* first_atom;
	Atom* temp_atom;
	Unit_Cell* next;
	double r[3];		//position vector of unit cell 

	Unit_Cell()		//Constructor
	{
		first_atom = NULL;
		temp_atom = NULL;
		next = NULL;
		r[0] = 0.0;
		r[1] = 0.0;
		r[2] = 0.0;
	}

	Unit_Cell(double r0_in, double r1_in, double r2_in)
	{
		first_atom = NULL;
		temp_atom = NULL;
		next = NULL;
		r[0] = r0_in;
		r[1] = r1_in;
		r[2] = r2_in;
	}


	void Add_Atom(Atom Atom_in)
	{
		if (first_atom == NULL)
		{
			first_atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.radius, Atom_in.m);
		}
		else
		{
			first_atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.radius, Atom_in.m);
		}
		//Give atoms absolute coordinates by adding displacement of current cell from the origin.

	}

	void Add_Atom(Atom* Atom_in)
	{
		if (first_atom == NULL)
		{
			first_atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m);
		}
		else
		{
			first_atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->radius, Atom_in->m);
		}

		//Give atoms absolute coordinates by adding displacement of current cell from the origin.

	}

	void Print()
	{
		temp_atom = first_atom;
		while (first_atom != NULL)
		{
			first_atom->Print();
			first_atom = first_atom->Next;
		}
		first_atom = temp_atom;
	}

	void Print(int lattice)
	{
		int atom=0;
		temp_atom = first_atom;
		while (first_atom != NULL)
		{
			first_atom->Print(lattice, atom);
			first_atom = first_atom->Next;
			atom++;
		}
		first_atom = temp_atom;
	}
};

class Lattice
{

public:

	Unit_Cell* first_cell;
	Unit_Cell* temp_cell;
	Unit_Cell* print_cell;
	Lattice* next;
	double a1[3];	//lattice vector 1
	double a2[3];	//lattice vector 2
	double a3[3];	//lattice vector 2
	double offset[3] = { 0,0,0 };

	Lattice(double a1_in[3], double a2_in[3], double a3_in[3])	//Constructor a1 and a2 are in nm, Phi is in degrees
	{
		first_cell = NULL;
		temp_cell = NULL;
		next = NULL;

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


	Lattice(double a1_in[3], double a2_in[3], double a3_in[3], double offset_in[3])	//Constructor a1 and a2 are in nm, Phi is in degrees
	{
		first_cell = NULL;
		temp_cell = NULL;
		next = NULL;

		a1[0] = a1_in[0];
		a1[1] = a1_in[1];
		a1[2] = a1_in[2];

		a2[0] = a2_in[0];
		a2[1] = a2_in[1];
		a2[2] = a2_in[2];

		a3[0] = a3_in[0];
		a3[1] = a3_in[1];
		a3[2] = a3_in[2];

		offset[0] = offset_in[0];
		offset[1] = offset_in[1];
		offset[2] = offset_in[2];
	}

	
	void Tile_Space(int a1_cells, int a2_cells, int a3_cells)
	{	
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
					first_cell= new Unit_Cell((a1[0] *(i+offset[0]) + a2[0] * (j + offset[1]) + a3[0] * (k + offset[2])), (a1[1] * (i + offset[0]) + a2[1] * (j + offset[1]) + a3[1] * (k + offset[2])) , (a1[2] * (i + offset[0]) + a2[2] * (j + offset[1]) + a3[2] * (k + offset[2])));
					first_cell->next = temp_cell;
				}
			}
		}
	}
	

	void Fill_Tiled_Space(Unit_Cell* Cell_in)
	{
		temp_cell = first_cell;
		Cell_in->temp_atom = Cell_in->first_atom;

		while (Cell_in->first_atom != NULL)
		{
			Cell_in->first_atom->r[0] = ((Cell_in->first_atom->r[0] * a1[0]) + (Cell_in->first_atom->r[1] * a2[0]) + (Cell_in->first_atom->r[2] * a3[0]));
			Cell_in->first_atom->r[1] = ((Cell_in->first_atom->r[0] * a1[1]) + (Cell_in->first_atom->r[1] * a2[1]) + (Cell_in->first_atom->r[2] * a3[1]));
			Cell_in->first_atom->r[2] = ((Cell_in->first_atom->r[0] * a1[2]) + (Cell_in->first_atom->r[1] * a2[2]) + (Cell_in->first_atom->r[2] * a3[2]));
			Cell_in->first_atom = Cell_in->first_atom->Next;
		}

		Cell_in->first_atom = Cell_in->temp_atom;

		while (first_cell->next != NULL)
		{
			while (Cell_in->first_atom != NULL)
			{		
				first_cell->Add_Atom(Cell_in->first_atom);
				Cell_in->first_atom = Cell_in->first_atom->Next;
			}
			Cell_in->first_atom = Cell_in->temp_atom;
			first_cell = first_cell->next;
		}
		first_cell = temp_cell;
	}

	void Print()
	{
		print_cell = first_cell;
		while (first_cell != NULL)
		{
			first_cell->Print();
			first_cell = first_cell->next;
		}
		first_cell = print_cell;
	}

	void Print(int lattice)
	{
		print_cell = first_cell;
		while (first_cell != NULL)
		{
			first_cell->Print(lattice);
			first_cell = first_cell->next;
		}
		first_cell = print_cell;
	}

	

};

class Surface
{
public:
	Lattice* first_Lattice;
	Lattice* temp_Lattice;
	Atom* first_Atom;
	Atom* temp_Atom;
	Atom* temp;
	int layer;

	Surface()		//Constructor
	{
		first_Lattice = NULL;
		first_Atom = NULL;
		temp = NULL;
	}

	Surface(Lattice* Lattice_in, Unit_Cell* Cell_in, int a1_cells, int	a2_cells, int a3_cells)
	{
		first_Lattice = Lattice_in;
		temp_Lattice = first_Lattice;
		while (first_Lattice != NULL)
		{
			first_Lattice->first_cell = Cell_in;
			first_Lattice->Tile_Space(a1_cells, a2_cells, a3_cells);
			first_Lattice->Fill_Tiled_Space(Cell_in);
			first_Lattice = first_Lattice->next;
		}
		first_Lattice = temp_Lattice;
		first_Atom = first_Lattice->first_cell->first_atom;
		temp = NULL;
	}


	void Print()
	{
		temp_Lattice = first_Lattice;
		while (first_Lattice != NULL)
		{
			first_Lattice->Print();
			first_Lattice = first_Lattice->next;
		}
		first_Lattice = temp_Lattice;
	}

	void Print_Lattice_Atom()
	{
		int lattice = 0;
		temp_Lattice = first_Lattice;
		while (first_Lattice != NULL)
		{
			first_Lattice->Print(lattice);
			first_Lattice = first_Lattice->next;
			lattice++;
		}
		first_Lattice = temp_Lattice;
	}

	
};

#endif