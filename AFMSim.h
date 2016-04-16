#ifndef AFMSim
#define AFMSim

const int LookupSize =128;
//Lookup table to be populated with Lennard-Jones parameters
double lookup[LookupSize][LookupSize][2] = {0};

int GlobalXRes = 500;
int GlobalYRes = 250;
int GlobalZRes = 250;

void write_xyz(string Filename)
{
	double x, y, z, dz;
    string   NewFile = Filename.substr(0,Filename.find_last_of(".")).append(".xyz");
    string   line;
    fstream  Input(Filename);
    ofstream Output(NewFile);
    stringstream StreamBuf;

    if (Input.is_open()) //skip header line
    {
        getline (Input,line);
        while ( Input >> x >> y >> z >> dz >>dz)
            StreamBuf << x/10 << ' ' << y/10 << ' ' << z/10 << '\n'; //Convert Angstroms to Nanometers

        Input.close();
    }
    else
        cout << "Unable to open file";

    Output << StreamBuf.str();
    Output.close();
}


class Atom
{
private:
    double temp;
    double e_mod;
    double s_mod;

    void populate_lookup()
    {
        if(lookup[Atomic_Number][Atomic_Number][0]!=0)
            return;

        else
        {
            lookup[Atomic_Number][Atomic_Number][0]=e;
            lookup[Atomic_Number][Atomic_Number][1]=s;

            for(int i=0; i<LookupSize; i++)
            {
                if(lookup[i][i][0]!=0 && i != Atomic_Number)
                {
                    lookup[Atomic_Number][i][0] =   sqrt(lookup[Atomic_Number][Atomic_Number][0]*lookup[i][i][0]);
                    lookup[Atomic_Number][i][1] =		(lookup[Atomic_Number][Atomic_Number][1]+lookup[i][i][1])/2;
                    lookup[i][Atomic_Number][0] =   	 lookup[Atomic_Number][i][0];
                    lookup[i][Atomic_Number][1] =   	 lookup[Atomic_Number][i][1];
                }
            }
        }
    }

public:
    double r[3];
    double	e;
    double	s;
    int Atomic_Number=0;
    string element;
    Atom* Next;

    Atom(double x_in, double y_in, double z_in, string element_in)
    {
        r[0] = x_in;
        r[1] = y_in;
        r[2] = z_in;
        element=element_in;

        if (element.compare("C") == 0)
        {
            Atomic_Number=6;
            e = 1.0456*pow(10, -21);		//epsilon	[Joules]	(LJ Parameter)
            s = 4.0;						//sigma		[Angstroms] (LJ Parameter)
        }
        else if (element.compare("Pt") == 0)
        {
            Atomic_Number=78;
            e = 8.3304*pow(10, -20);		//epsilon	[Joules]	(LJ Parameter)
            s = 2.475;						//sigma		[Angstroms] (LJ Parameter)
        }

        else if (element.compare("N") == 0)
        {
            Atomic_Number=7;
            e = 8.1111631*pow(10, -21);		//epsilon	[Joules]	(LJ Parameter)
            s = 3.5;						//sigma		[Angstroms] (LJ Parameter)
        }

        else if (element.compare("O") == 0)
        {
            Atomic_Number=8;
            e = 1.38953*pow(10, -21);		//epsilon	[Joules]	(LJ Parameter)
            s = 3.20;						//sigma		[Angstroms] (LJ Parameter)
        }

        else if (element.compare("S") == 0)
        {
            Atomic_Number=16;
            e = 1.38953*pow(10, -21);		//epsilon	[Joules]	(LJ Parameter)
            s = 4.00;						//sigma		[Angstroms] (LJ Parameter)
        }

        else if (element.compare("H") == 0)
        {
            Atomic_Number=1;
            e = 1.38953*pow(10, -22);		//epsilon	[Joules]	(LJ Parameter)
            s = 2.00;						//sigma		[Angstroms] (LJ Parameter)
        }


        else
        {
            cout << "Lennard jones parameters not stored for: " << element << ". Assuming Carbon"<< endl;
            element = "C";
            Atomic_Number=8;
            e = 1.0456*pow(10, -21);	//epsilon	[Joules]	(LJ Parameter)
            s = 4.0;					//sigma		[Angstroms] (LJ Parameter)
        }

        populate_lookup();
        Next = NULL;
    }


    void Print()
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
    }

    void Print(int var)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var << "\n";
    }

    void Print(double var)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var << "\n";
    }

    void Print(double var1, double var2)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var1 << "\t" << var2 << "\n";
    }

    double Dist(Atom* Other_Atom)
    {
        return (double)  sqrt(pow((Other_Atom->r[2] - r[2]), 2) + pow((Other_Atom->r[1] - r[1]), 2) + pow((Other_Atom->r[0] - r[0]), 2));
    }

    double LJPot(Atom* Other_Atom) // Calculate interaptomic potemtial
    {
        temp = Dist(Other_Atom);
        if (temp > (3 * lookup[Atomic_Number][Other_Atom->Atomic_Number][1])) 	// Truncation outside critical range
            return 0;

        e_mod= lookup[Atomic_Number][Other_Atom->Atomic_Number][0];
        s_mod= lookup[Atomic_Number][Other_Atom->Atomic_Number][1];

        const double LJ_Truncation = 4 * e_mod *(pow((s_mod / (3 * s_mod)), 12) - pow((s_mod / (3 * s_mod)), 6));

        return (double)	 4 * e_mod *(pow((s_mod /temp), 12) - pow((s_mod / temp), 6)) - LJ_Truncation;
    }

    double LJForce(Atom* Other_Atom) // Calculate interatomic force
    {
        temp = Dist(Other_Atom);
        if (temp > (3 * lookup[Atomic_Number][Other_Atom->Atomic_Number][1])) 	// Truncation outside critical range
            return 0;

        e_mod= lookup[Atomic_Number][Other_Atom->Atomic_Number][0];
        s_mod= lookup[Atomic_Number][Other_Atom->Atomic_Number][1];

        const double LJ_Truncation = (24 / (3*s_mod))*e_mod*(2 * pow((s_mod / (3 * s_mod)), 12)- pow((s_mod / (3 * s_mod)), 6))*pow(10, 19);

        return (double) (24 / temp)*e_mod*(2 * pow((s_mod / temp), 12) -  pow((s_mod / temp), 6))*pow(10, 19); //- LJ_Truncation;
    }

    double PerpLJForce(Atom* Other_Atom) // Calculate interatomic force lying in z direction
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
            atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.element);

        else
        {
            while (atom->Next != NULL)
            {
                atom = atom->Next;
            }
            atom_count++;
            atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.element);
            atom = Temp;
        }
    }

    void Add_Atom(Atom* Atom_in)
    {
        Atom* Temp = atom;
        if (atom == NULL)
            atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->element);

        else
        {
            while (atom->Next != NULL)
            {
                atom = atom->Next;
            }
            atom_count++;
            atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2],Atom_in->element);
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

    void Import(string filename_in) //Import vesta object from .vesta AND .xyz file
    {
        string line;
        string vesta_filename = filename_in;
        string xyz_filename = filename_in;
        string element;
        ifstream fin;
        int i = 1;
        double junk;
        double x,y,z;
        double initial[3];
        double rotation[3][3];
        double rotation2[3][3] =
        {
            { 1,	0,	 0},
            { 0,	0, 	-1},
            { 0,	1,	 0},
        };


        Atom* Temp;

        ResetOrigin();

        for (int i = 0; i < 3; i++)
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

        for (int i = 0; i < 3; i++)
        {
            fin >> rotation[i][0]  >> rotation[i][1]  >> rotation[i][2]  >> junk;
        }

        fin.close();

        fin.open(xyz_filename);
        i=0;
        if (fin.is_open())
        {
            getline(fin, line);
            getline(fin, line);

            while (fin >> element >> x >> y >> z)
            {
                i++;
                Temp = new Atom(x,y,z,element);
                cout << "Atom " << setw(2) << i << ": " << Temp->element << "\tZ= " <<Temp->Atomic_Number<<endl;

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

        for (int i = 0; i < 3; i++)
            r[i] = initial[i];

        ResetOrigin();

    }

    void RotateAboutX(double Ang)
    {
        Ang= M_PI*Ang/180;
        double Roatation[3][3] =
        {
            { 1,		0,				0		},
            { 0,	cos(Ang),	-1.0 *sin(Ang)	},
            { 0,	sin(Ang),		  cos(Ang)	},

        };
        Rotate(Roatation);
    }

    void RotateAboutY(double Ang)
    {
        Ang= M_PI*Ang/180;
        double Roatation[3][3] =
        {
            { cos(Ang),			0,		sin(Ang) },
            { 0,				1,		0		 },
            { sin(Ang)*-1.0,	0,		cos(Ang) },
        };
        Rotate(Roatation);
    }

    void RotateAboutZ(double Ang)
    {
        Ang= M_PI*Ang/180;
        double Roatation[3][3] =
        {
            { cos(Ang),	-1.0*	sin(Ang),	0 },
            { sin(Ang),			cos(Ang),	0 },
            { 0,				0,			1 },
        };
        Rotate(Roatation);
    }

    void Rotate(double rotation[3][3])
    {
        double r_rot[3];
        Atom* temp = atom;
        double inintial[3] = { r[0],r[1],r[2] };
        Move(-1.0 * r[0], -1.0 * r[1], -1.0 * r[2]);

        while (temp != NULL)
        {
            r_rot[0] = rotation[0][0] * temp->r[0] + rotation[0][1] * temp->r[1] + rotation[0][2] * temp->r[2];
            r_rot[1] = rotation[1][0] * temp->r[0] + rotation[1][1] * temp->r[1] + rotation[1][2] * temp->r[2];
            r_rot[2] = rotation[2][0] * temp->r[0] + rotation[2][1] * temp->r[1] + rotation[2][2] * temp->r[2];

            temp->r[0] = r_rot[0];
            temp->r[1] = r_rot[1];
            temp->r[2] = r_rot[2];

            temp = temp->Next;
        }
        Move(inintial[0], inintial[1], inintial[2]);
    }

    void ResetOrigin() //Set lowest point on tip as coordinate
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
        r[2] = 0;

        atom = temp;
    }

    void Move(double x_in, double y_in, double z_in) //Move tip by x,y,z
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
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\n";
    }

    void Print(int var)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var << "\n";
    }

    void Print(double var)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var << "\n";
    }

    void Print(double var1, double var2)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var1 << "\t" << var2 << "\n";
    }

    void Print(double var1, double var2, double var3)
    {
        cout << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << var1 << "\t" << var2 << "\t" << var3 << "\n";
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
            first_atom = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.element);

        else
        {
            while (first_atom->Next != NULL)
            {
                first_atom = first_atom->Next;
            }
            first_atom->Next = new Atom(Atom_in.r[0] + r[0], Atom_in.r[1] + r[1], Atom_in.r[2] + r[2], Atom_in.element);
            first_atom = Temp;
        }
    }

    void Add_Atom(Atom* Atom_in)
    {
        Atom* Temp = first_atom;
        if (first_atom == NULL)
            first_atom = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->element);

        else
        {
            while (first_atom->Next != NULL)
            {
                first_atom = first_atom->Next;
            }
            first_atom->Next = new Atom(Atom_in->r[0] + r[0], Atom_in->r[1] + r[1], Atom_in->r[2] + r[2], Atom_in->element);
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


    void Tile_Space() //Fill space with empty unit cells until surface of desired size is generated
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

    void Fill_Tiled_Space(Unit_Cell* Cell_in) //Fill empty unit cells with atoms as in Cell_in
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

    //Control loop ajusting tip height until desired force is found, similar to real AFM operation
    void TipHeightMap(VestaObject* Obj_in, double XMax, double XMin, double YMax, double YMin, double ZIni, double Setpoint, double ZStep, double FRes)
    {
        int progress = 0;
        double Prev = 0;
        double Cur = 0;
        double xstep = (XMax - XMin) / GlobalXRes;
        double ystep = (YMax - YMin) / GlobalYRes;

        if (XMax == XMin)
            xstep = 1;

        if (YMax == YMin)
            ystep = 1;

        double original_position[3] = { Obj_in->r[0], Obj_in->r[1], Obj_in->r[2] };

        cout << "x[Angstroms]\t" << "y[Angstroms]\t" << "z[Angstroms]\t" << "dz\t" << endl;

        Obj_in->Move(XMin - Obj_in->r[0], YMin - Obj_in->r[1], ZIni - Obj_in->r[2]);


        while (Obj_in->r[0] <= XMax)
        {

            Prev = TipHeightCalc(Obj_in, Setpoint, ZStep, FRes);
            Cur = TipHeightCalc(Obj_in, Setpoint, ZStep, FRes);
            while (Obj_in->r[1] <= YMax)
            {
                Cur  = TipHeightCalc(Obj_in, Setpoint, ZStep, FRes);
                Obj_in->Print((Cur-Prev),Setpoint);
                Prev = Cur;
                Obj_in->Move(0, ystep, 0);
            }
            Obj_in->Move(xstep, YMin - Obj_in->r[1], 0);

            progress = (int)((Obj_in->r[0] - XMin) * 100 / (XMax - XMin));
            printf("\rTopology %d%% Complete", progress);
            fflush( stdout );
            if (YMax != YMin)
            {
                cout << "\n";
            }
        }
        printf("\n");
        Obj_in->Move(original_position[0] - Obj_in->r[0], original_position[1] - Obj_in->r[1], original_position[2] - Obj_in->r[2]);
    }

    //Function responsible for chosing direction of tip movement in control loop
    double TipHeightCalc(VestaObject* Obj_in, double Setpoint, double ZStep, double FRes)
    {
        double Error = 100000;
        bool check1 = 0;
        bool check2 = 0;

        while (abs(LJForce(Obj_in)) < 0.001)
            Obj_in->Move(0, 0, -1.0*ZStep);


        while (abs(Error) >  FRes)
        {
            while ((check1 * check2) != 1)
            {
                Error = LJForce(Obj_in) - Setpoint;

                if (Error > 0)
                {
                    Obj_in->Move(0, 0,ZStep);
                    check1 = 1;
                }

                else
                {
                    Obj_in->Move(0,0, -1.0*ZStep);
                    check2 = 1;
                }
            }

            ZStep = ZStep /5;
            check1 = 0;
            check2 = 0;
        }
        return(Obj_in->r[2]);
    }

public:

    int cell_count = 1;
    int atom_count = 0;
    double a1[3];		//////////////////////////////////////////////////////////////////
    double a2[3];		//		Lattice vectors a1, a2 and a3 follow from
    double a3[3];		//	crystallographic convention and have units of Angstroms
    int a1_cells;		//////////////////////////////////////////////////////////////////
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

    double LJPot(VestaObject* Obj_in) //Calculate potential between tip and surface
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

    double LJForce(VestaObject* Obj_in) //Calculate Force between tip and surface
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

    double PerpLJForce(VestaObject* Obj_in) //Calculate force between tip and surface in z direction
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

    //Select area to be scanned for topology measurement
    void TipHeight(VestaObject* Obj_in, double Setpoint, double ZRes, double FRes, double Area)
    {
        double XMax = round(((a1[0] * a1_cells + a2[0] * a2_cells) / 2 + Area * (abs(a1[0]) + abs(a2[0]))));
        double XMin = round(((a1[0] * a1_cells + a2[0] * a2_cells) / 2 - Area * (abs(a1[0]) + abs(a2[0]))));
        double YMax = round(((a1[1] * a1_cells + a2[1] * a2_cells) / 2 + Area * (abs(a1[1]) + abs(a2[1]))));
        double YMin = round(((a1[1] * a1_cells + a2[1] * a2_cells) / 2 - Area * (abs(a1[1]) + abs(a2[1]))));
        std::stringstream buffer;
        std::ofstream out("Topology.dat");

        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        std::cout.rdbuf(buffer.rdbuf());

        TipHeightMap(Obj_in, XMax, XMin, YMax, YMin, Obj_in->r[2], Setpoint, ZRes, FRes);

        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat
        cout << buffer.str() << endl;
        std::cout.rdbuf(coutbuf); //reset to standard output again
    }

    //Select area to be scanned for XY Plane force measurement
    void SurfaceForce(VestaObject* Obj_in, double z_0, double Area)
    {
        double XMax = round((a1[0] * a1_cells + a2[0] * a2_cells) / 2 + Area * (abs(a1[0]) + abs(a2[0])));
        double XMin = round((a1[0] * a1_cells + a2[0] * a2_cells) / 2 - Area * (abs(a1[0]) + abs(a2[0])));
        double YMax = round((a1[1] * a1_cells + a2[1] * a2_cells) / 2 + Area * (abs(a1[1]) + abs(a2[1])));
        double YMin = round((a1[1] * a1_cells + a2[1] * a2_cells) / 2 - Area * (abs(a1[1]) + abs(a2[1])));
        std::stringstream buffer;
        std::ofstream out("Surface.dat");

        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        std::cout.rdbuf(buffer.rdbuf());

        XYZForce(Obj_in, XMax, XMin, YMax, YMin, z_0, z_0, 1);

        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to Surface.dat
        cout << buffer.str() << endl;
        std::cout.rdbuf(coutbuf); //reset to standard output again
    }

    //Select area to be scanned for XZ Plane force measurement
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

    //Function used to scan tip over a volume, called by Forcecurve() and SurfaceForce()
    void XYZForce(VestaObject* Obj_in, double XMax, double XMin, double YMax, double YMin, double ZMax, double ZMin, bool Corrugation)
    {
        double Average = 0;
        int progress = 0;
        int i = 0;
        double xstep = (XMax - XMin) / GlobalXRes;
        double ystep = (YMax - YMin) / GlobalYRes;
        double zstep = (ZMax - ZMin) / GlobalZRes;

        if (XMax == XMin)
            xstep = 1;

        if (YMax == YMin)
            ystep = 1;

        if (ZMax == ZMin)
            zstep = 1;

        double original_position[3] = { Obj_in->r[0], Obj_in->r[1], Obj_in->r[2] };

        cout << "x[Angstroms]\t" << "y[Angstroms]\t" << "z[Angstroms]\t" << "F[nN]\t" << endl;

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
                    printf("\rForce calculation %d%% Complete", progress);
                    fflush( stdout );
                }
            }

            if (ZMax > ZMin*1.01)
            {
                progress = (int)((Obj_in->r[2] - ZMin) * 100 / (ZMax - ZMin));
                printf("\rForce calculation %d%% Complete", progress);
                fflush( stdout );

            }

            Obj_in->Move(XMin- Obj_in->r[0], YMin - Obj_in->r[1], zstep);
            cout << "\n";
            Average = 0;
            i = 0;
        }
        printf("\n");
        Obj_in->Move(original_position[0] - Obj_in->r[0], original_position[1] - Obj_in->r[1], original_position[2] - Obj_in->r[2]);
    }

    void Add_Defect(Unit_Cell* Defect_Cell) //Add vesta object to surface
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

    void Remove_Defect()
    {
        Unit_Cell* Temp = first_cell;
        while (first_cell->Next->Next != NULL)
            first_cell = first_cell->Next;

        first_cell->Next = NULL;
        first_cell = Temp;
    }
};

//Work in progress, not yet correct beyond this point
class Cantilever
{
private:
    typedef boost::array< double, 2 > state_type;
    typedef runge_kutta_cash_karp54< state_type > stepper_type;
    double ft[2048];
    double z_prev = 0;

public:
    double l;
    double w;
    double t;
    double k;
    double c = 0;
    double Fd = 0; //driving force;
    double rho;
    double E;
    double I;
    double m;
    double w_0;
    double f_0;
    double Z_0;
    int i = 0;
    const int transient = 1000;
    VestaObject* Tip;
    Surface* Surf;

    Cantilever(double length, double width, double thickness, double Youngs_Modulus, double density, VestaObject* Tip_in, Surface* Surf_in)
    {
        l = length;
        w = width;
        t = thickness;
        E = Youngs_Modulus;
        I = l*pow(t, 3) / 12;
        k = 3 * E*I / pow(l, 3);
        rho = density;
        m = l*w*t*rho;
        w_0 = sqrt(k / m);
        f_0 = w_0 / (2 * M_PI);
        Tip = Tip_in;
        Surf = Surf_in;
        cout << "f_0: " << f_0 << endl;
        cout << "k: " << k << endl;

    }

    void four1(double* data, unsigned long nn)
    {
        unsigned long n, mmax, m, j, istep, i;
        double wtemp, wr, wpr, wpi, wi, theta;
        double tempr, tempi;

        // reverse-binary reindexing
        n = nn << 1;
        j = 1;
        for (i = 1; i<n; i += 2)
        {
            if (j>i)
            {
                swap(data[j - 1], data[i - 1]);
                swap(data[j], data[i]);
            }
            m = nn;
            while (m >= 2 && j>m)
            {
                j -= m;
                m >>= 1;
            }
            j += m;
        };

        // here begins the Danielson-Lanczos section
        mmax = 2;
        while (n>mmax)
        {
            istep = mmax << 1;
            theta = -(2 * M_PI / mmax);
            wtemp = sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            for (m = 1; m < mmax; m += 2)
            {
                for (i = m; i <= n; i += istep)
                {
                    j = i + mmax;
                    tempr = wr*data[j - 1] - wi*data[j];
                    tempi = wr * data[j] + wi*data[j - 1];

                    data[j - 1] = data[i - 1] - tempr;
                    data[j] = data[i] - tempi;
                    data[i - 1] += tempr;
                    data[i] += tempi;
                }
                wtemp = wr;
                wr += wr*wpr - wi*wpi;
                wi += wi*wpr + wtemp*wpi;
            }
            mmax = istep;
        }
    }

    void Oscillate(const state_type &z, state_type &dzdt, double t)
    {
        dzdt[0] = z[1];
        dzdt[1] = ((-1.0*k)*z[0] - c*dzdt[0] + Fd*sin(w_0*t)  + Surf->LJForce(Tip)*pow(10, -9))/m; //restoring force - damping + driving + Tip-Surface interaction
    }

    void Write_Oscillation(const state_type &z, const double t)
    {

        Tip->Print(t, z[0], z[1]);
        Tip->Move(0, 0, (z[0]*pow(10,10) - z_prev*pow(10, 10)));
        z_prev = z[0];
        if (t > transient / f_0 && i < 1024)
        {
            ft[i] = z[0];
            i++;
        }
    }

    void Write_Bifurcation(const state_type &z, const double t)
    {
        if (t > (transient + i)/f_0)
        {
            Tip->Print(t, z[0], Z_0);
            i++;
            Tip->Move(0, 0, (z[0] - z_prev)*pow(10, 10));
            z_prev = z[0];
        }
    }

    void FM_Scan(Surface* Surf, double Area)
    {
        double XMax = round(((Surf->a1[0] * Surf->a1_cells + Surf->a2[0] * Surf->a2_cells) / 2 + Area * (abs(Surf->a1[0]) + abs(Surf->a2[0]))));
        double XMin = round(((Surf->a1[0] * Surf->a1_cells + Surf->a2[0] * Surf->a2_cells) / 2 - Area * (abs(Surf->a1[0]) + abs(Surf->a2[0]))));
        double YMax = round(((Surf->a1[1] * Surf->a1_cells + Surf->a2[1] * Surf->a2_cells) / 2 + Area * (abs(Surf->a1[1]) + abs(Surf->a2[1]))));
        double YMin = round(((Surf->a1[1] * Surf->a1_cells + Surf->a2[1] * Surf->a2_cells) / 2 - Area * (abs(Surf->a1[1]) + abs(Surf->a2[1]))));
        double XStep = (XMax - XMin) / 250;
        double YStep = (YMax - YMin) / 250;
        double Original_Position[3] = { Tip->r[0], Tip->r[1], Tip->r[2] };
        namespace pl = std::placeholders;
        const double end = (transient+50)/f_0;
        const double step = 0.00005 / f_0;


        Tip->Move(XMin - Tip->r[0], YMin - Tip->r[1], 0);

        std::stringstream buffer;
        std::ofstream Osc("Oscillation.dat");
        std::ofstream Transform("Transform.dat");
        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        std::cout.rdbuf(Osc.rdbuf()); //redirect std::cout to Oscillation.dat

        state_type z = { 0 , 0 }; // initial conditions
        c  = 1000*m;
        Fd = m/100;

        integrate_adaptive(make_controlled(1E-12 , 1E-12 ,stepper_type()),
                           std::bind(&Cantilever::Oscillate, *this, pl::_1, pl::_2, pl::_3), z, 0.0, end, step, std::bind(&Cantilever::Write_Oscillation, *this, pl::_1, pl::_2));


        std::cout.rdbuf(Transform.rdbuf()); //redirect std::cout to Transform.dat"


        four1(ft, 1024);
        double conv = 2*1024 * step;
        for (int k = 0; k < 2047; k+=2)
        {
            cout << k / conv << " " << sqrt(pow(ft[k], 2) + pow(ft[k + 1], 2)) << endl;
        }

        std::cout.rdbuf(coutbuf); //reset to standard output again
    }

};

#endif