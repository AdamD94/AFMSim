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
	double r; 
	double m;
	Atom* Next;

	Atom(float x_in, float y_in, float r_in, float m_in)
	{
		x = x_in;
		y = y_in;
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

		Lattice(float mag_a1_in, float mag_a2_in, float Phi_in)	//Constructor
		{
			
			Phi = Phi_in * (M_PI/180.0);

			a1[0] = mag_a1_in * pow(10, -9);	//define a1 as lying on the x-axis;
			a1[1] = 0;							// Zero Y component;
			a1[2] = 0;							// Zero Z component;

			a2[0] = mag_a2_in * pow(10, -9) * cos(180.0 - Phi);	//x coordinate of a2
			a2[1] = mag_a2_in * pow(10, -9) * sin(180.0 - Phi);	//y coordinate of a2
			a2[2] = 0;											//z coordinate of a2
		}
	};
	
	class Basis
	{
	public:
		Atom* first_atom;
			

		Basis()		//Constructor
		{
			first_atom = NULL;
		}
	};          

public:
	Basis*	first_basis;



	Surface()		//Constructor
	{
		first_basis = NULL;
	}

};

int	main()
{


	return(0);
}
