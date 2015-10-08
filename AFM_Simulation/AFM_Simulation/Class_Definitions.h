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

};

class Unit_Cell
{
public:
	class Lattice;
	Atom* first_atom;
	Atom* final_atom;
	Atom* temp_atom;
	Unit_Cell* next;
	double r[3];		//position vector of unit cell 

	Unit_Cell()		//Constructor
	{
		first_atom = NULL;
		final_atom = NULL;
		temp_atom = NULL;
		next = NULL;
		r[0] = 0.0;
		r[1] = 0.0;
		r[2] = 0.0;
	}

	Unit_Cell(double r0_in, double r1_in, double r2_in)
	{
		first_atom = NULL;
		final_atom = NULL;
		temp_atom = NULL;
		next = NULL;
		r[0] = r0_in;
		r[1] = r1_in;
		r[2] = r2_in;
	}


	void Add_Atom(Atom Atom_in)
	{
		//temp_atom = first_atom;
		//first_atom = new Atom(Atom_in.r[0], Atom_in.r[1], Atom_in.r[2], Atom_in.radius, Atom_in.m);
		//first_atom->Next = temp_atom;
		if (first_atom == NULL)
		{
			first_atom = new Atom(Atom_in.r[0], Atom_in.r[1], Atom_in.r[2], Atom_in.radius, Atom_in.m);
		}
		else
		{
			first_atom->Next = new Atom(Atom_in.r[0], Atom_in.r[1], Atom_in.r[2], Atom_in.radius, Atom_in.m);
		}
		//Give atoms absolute coordinates by adding displacement of current cell from the origin.
		first_atom->r[0] += this->r[0];
		first_atom->r[1] += this->r[1];
		first_atom->r[2] += this->r[2];
		cout << first_atom->r[0] << " " << first_atom->r[1] << " " << first_atom->r[2] <<  endl;
	}

	void Add_Atom(Atom* Atom_in)
	{
		//temp_atom = first_atom;
		//first_atom = new Atom(Atom_in->r[0], Atom_in->r[1], Atom_in->r[2], Atom_in->radius, Atom_in->m);
		//first_atom->Next = temp_atom;
		if (first_atom == NULL)
		{
			first_atom = new Atom(Atom_in->r[0], Atom_in->r[1], Atom_in->r[2], Atom_in->radius, Atom_in->m);
		}
		else
		{
			first_atom->Next = new Atom(Atom_in->r[0], Atom_in->r[1], Atom_in->r[2], Atom_in->radius, Atom_in->m);
		}

		//Give atoms absolute coordinates by adding displacement of current cell from the origin.
		first_atom->r[0] += this->r[0];
		first_atom->r[1] += this->r[1];
		first_atom->r[2] += this->r[2];
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

	


};

class Lattice
{

public:

	Unit_Cell* first_cell;
	Unit_Cell* final_cell;
	Unit_Cell* temp_cell;
	Unit_Cell* temp_cell_2;
	Lattice* next;
	double a1[3];	//lattice vector 1
	double a2[3];	//lattice vector 2

	Lattice(double a1_in[3], double a2_in[3])	//Constructor a1 and a2 are in nm, Phi is in degrees
	{
		first_cell = NULL;
		final_cell = NULL;
		temp_cell = NULL;
		temp_cell_2 = NULL;
		next = NULL;

		a1[0] = a1_in[0];
		a1[1] = a1_in[1];
		a1[2] = a1_in[2];

		a2[0] = a2_in[0];
		a2[1] = a2_in[1];
		a2[2] = a2_in[2];
	}

	void Generate_Lattice(double a1_in[3], double a2_in[3])
	{
		first_cell = NULL;
		final_cell = NULL;
		temp_cell = NULL;
		next = NULL;

		a1[0] = a1_in[0];
		a1[1] = a1_in[1];
		a1[2] = a1_in[2];

		a2[0] = a2_in[0];
		a2[1] = a2_in[1];
		a2[2] = a2_in[2];
	}
	
	void Tile_Space(double x_max, double y_max)
	{
		while (first_cell->r[1] <= (y_max))
		{
			while (first_cell->r[0] <= (x_max ))
			{
				temp_cell = first_cell;
				first_cell = new Unit_Cell((temp_cell->r[0] + a1[0]), (temp_cell->r[1] + a1[1]), (temp_cell->r[2] + a1[2]));
				first_cell->next = temp_cell;
			}
			first_cell->r[0] = 0;
			temp_cell = first_cell;
			first_cell = new Unit_Cell((temp_cell->r[0] + a2[0]), (temp_cell->r[1] + a2[1]), (temp_cell->r[2] + a2[2]));
			first_cell->next = temp_cell;
		}
	}
	

	void Fill_Tiled_Space(Unit_Cell* Cell_in)
	{
		temp_cell = first_cell;
		Cell_in->temp_atom = Cell_in->first_atom;
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
		temp_cell = first_cell;
		while (first_cell != NULL)
		{
			first_cell->Print();
			first_cell = first_cell->next;
		}
		first_cell = temp_cell;
	}
};

class Surface
{
public:
	Lattice* first_Lattice;
	Atom* first_Atom;
	Atom* temp_Atom;
	Atom* temp;

	Surface()		//Constructor
	{
		first_Lattice = NULL;
		first_Atom = NULL;
		temp = NULL;
	}

	void Generate_Surface(Lattice* Lattice_in, Unit_Cell* Cell_in, double x_max, double y_max)
	{
		first_Lattice = Lattice_in;
		first_Lattice->first_cell = Cell_in;
		first_Lattice->Tile_Space(x_max, y_max);
		first_Lattice->Fill_Tiled_Space(Cell_in);
		first_Atom = first_Lattice->first_cell->first_atom;
	}

	void Print()
	{
	
		first_Lattice->Print();

	}

};

#endif