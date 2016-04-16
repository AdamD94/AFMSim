#pragma once
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <bitset>
#include "boost/filesystem.hpp"
#include "boost/array.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/range.hpp"
#include "gnuplot-iostream.h"
using namespace boost::numeric::odeint;
using namespace boost::filesystem;
using namespace std;

bool Defect_Present = 0;

int a1_cells 		= 30;		// Number of cells to be generated in the a1, a2 and a3 directions
int a2_cells 		= 30;
int a3_cells 		= 1;

double a 			= 1.4;
double c 			= 6.71;
double z_0			= 4;
double ZStep		= 0.001;
double FRes			= 0.001;	//Tolerance on force measurement
double Setpoint 	= 0;		//Force to be followed by tip nanoNewtons
double Ang			= 0;
double a1[3] 		= { a* 3.0 / 2.0	, a * sqrt(3.0) / 2.0	,a*  0.0 };	// Lattice vectors as in crystallography
double a2[3] 		= { a* 3.0 / 2.0	, a *-sqrt(3.0) / 2.0	,a*  0.0 };
double a3[3]		= { a				, 0						,c* -1.0 };

string Defect_Filename 	="";
string Tip_Filename 	="";
string Functionalise	="";