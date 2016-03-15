#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <direct.h>
#include "boost/array.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/range.hpp"

using namespace std;
using namespace boost::numeric::odeint;
typedef boost::array< double, 2 > state_type;
typedef runge_kutta_cash_karp54< state_type > stepper_type;

#include "Class_Definitions.h"  //Header File containing class definitions

double a = 1.4;
double c = 6.71;
double z_0 = 4;
double ZStep = 0.001;
double FRes = 0.00005; //Tolerance on force measurement
double Setpoint = 0; //Force to be followed by tip nanoNewtons
double Ang = 0;
int a1_cells = 40;		// Number of cells to be generated in the a1, a2 and a3 directions
int	a2_cells = 40;
int a3_cells = 1;

double a1[3] = { a* 3.0 / 2.0	, a * sqrt(3.0) / 2.0	,a*  0.0 };	// Lattice vectors as in crystallagraphy
double a2[3] = { a* 3.0 / 2.0	, a *-sqrt(3.0) / 2.0	,a*  0.0 };
double a3[3] = { a				, 0						,c* -1.0 };
bool Defect_Present = 0;
string Defect_Filename = "";
string Tip_Filename = "";
VestaObject* Defect = new VestaObject(0, 0, 0);
