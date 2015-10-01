using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

class Atom
{
public:
	double x;
	double y;
	double z;
	double r; 
	double m;
	Atom* Next;

	Atom(float x_in, float y_in, float z_in, float r_in, float m_in)
	{
		x = x_in;
		y = y_in;
		z = z_in;
		r = r_in;
		m = m_in;
		Next = NULL;
	}
};

class Surface
{
	class Lattice
	{
	public:
		double a1[3];	//lattice vector 1
		double a2[3];	//lattice vector 2
		double Phi;		//Angle between lattice vectors in radians

		Lattice(double mag_a1_in, double mag_a2_in, double Phi_in, double z_displacement)	//Constructor a1 and a2 are in nm, Phi is in degrees
		{
			
			Phi = Phi_in * (M_PI/180.0);

			a1[0] = mag_a1_in * pow(10, -9);		// Define a1 as lying on the x-axis;
			a1[1] = 0;								// Zero Y component by definition;
			a1[2] = z_displacement * pow(10, -9);	// Height above z=0;

			a2[0] = mag_a2_in * pow(10, -9) * cos(180.0 - Phi);	//x coordinate of a2
			a2[1] = mag_a2_in * pow(10, -9) * sin(180.0 - Phi);	//y coordinate of a2
			a2[2] = z_displacement * pow(10, -9);				//z coordinate of a2
		}
	};
	
	class Basis
	{
	public:
		Atom* first_atom;
		Atom* final_atom;
		Basis* next_cell;
		Basis* prev_cell;
		double x;
		double y;
		double z;

		Basis()		//Constructor
		{
			first_atom = NULL;
			final_atom = NULL;
			next_cell = NULL;
			prev_cell = NULL;
			x = 0;
			y = 0;
			z = 0;
		}

		void Add_Atom(Atom Atom_in)
		{			
			if (first_atom = NULL)
			{
				first_atom = &Atom_in;
				first_atom->x += this->x;
				first_atom->y += this->y;
				first_atom->z += this->z;

				final_atom = &Atom_in;
			}
			else
			{
				final_atom->Next = &Atom_in;
				final_atom = final_atom->Next;
				final_atom->x += this->x;
				final_atom->y += this->y;
				final_atom->z += this->z;
			}

			
		}
	};          

public:
	Basis*	first_basis;
	Lattice* first_lattice;

	Surface()		//Constructor
	{
		first_basis = NULL;
		first_lattice = NULL;
	}

};

int	main()
{

	Atom Carbon(0, 0, 0, 1.0, 1.0);  //x,y,z,r,m;
	Surface* Square = new Surface();
	Square->first_basis->Add_Atom(Carbon);
	return(0);
}
