#include <iostream>
#include<fstream>
#include<cmath>
#include<memory>
#include<stdlib.h>
#include <sys/time.h>
#include "Grid.h"
#include "FileReader.h"
#include "VTKFileWriter.h"

using namespace std;

typedef double Real;

#define D 2
#define Q 9
//Global Variables
size_t sizex, sizey, timesteps;
Real omega;
bool flag = false;
Grid *fluid, *tmpfluid;
double uw[2] = {0.08, 0.0};
double disvel[Q][2] = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 }, { -1.0, 1.0 }, { -1.0, 0.0 }, { -1.0, -1.0 }, { 0.0, -1.0 }, { 1.0, -1.0 } };
int neighbours[Q][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 },  { 0, -1 }, { 1, -1 } };
double stencil[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 9.0, 1.0 / 36.0 };
string imfile;

inline void init()
{
	ifstream infile;
	size_t * arr;
	if(flag)
	{
	 infile.open(imfile);
	 stringstream ss;
	 string inputLine = "";
        
        // Check if pgm-version is P2 (graymap)
        getline(infile,inputLine);
        if(inputLine.compare("P2") != 0) cerr << "Error: PGM-File is not a graymap" << endl;

        // Ignore second line (comment) in pgm-file
        getline(infile,inputLine);

        // Continue with a stringstream
        ss << infile.rdbuf();
        // Get size of domain
        ss >> sizex >> sizey;
        cout << sizex << " columns and " << sizey << " rows" << endl;
        
        // allocate array to store pixel data
        arr = new size_t[sizex * sizey];
        size_t value;
        ss >> value;
        
        // Iterate through pixel data
        for(size_t it_col = 0; it_col < sizey; ++it_col)
        {
                for (size_t it_row = 0; it_row < sizex ; ++it_row)
                {
                        ss >> arr[it_col * sizex + it_row];
                        //cout << arr[it_col * sizex +  it_row] << " ";
                }
                //cout << endl;
        }     
        // close file
        infile.close();
	}
	double ux = 0.0, uy = 0.0, initrho =1.0;
	fluid = new Grid(sizex+2, sizey+2, stencil, ux, uy, initrho, arr , flag);
	tmpfluid = new Grid(sizex+2, sizey+2, stencil, ux, uy, initrho, arr , flag);

}

inline double calLidVel(size_t k)
{
	return (6.0 * stencil[k] * (disvel[k][0] * uw[0] + disvel[k][1] * uw[1]));
}

inline void stream()
{
	for (int i = 1; i <= (int)sizey; i++)
	{
		for (int j = 1; j <= (int)sizex; j++)
		{
        if(flag && (*tmpfluid).getBoundary(i, j) == 0)
		{	
		for (int k = 1; k < Q; k++)
			{
                int pt = (*tmpfluid).getBoundary(i + neighbours[k][1], j + neighbours[k][0]);
				size_t l;
                if(k != 4) l = (k+(Q-1)/2)%(Q-1);
                else l=8;

                if (pt > 0 && pt != 3 )
			(*tmpfluid)(i, j, l) = (*fluid)(i, j, k);       // Bounce back from vertical and bottom walls					
                else if (pt == 3)
			(*tmpfluid)(i, j, l) = (*fluid)(i, j, k) - calLidVel(k); // Streaming effect due to lid velocity in upper wall
		else
                    (*tmpfluid)(i + neighbours[k][1], j + neighbours[k][0], k) = (*fluid)(i,j,k);  // Free Streaming
				
			}

			(*tmpfluid)(i, j, 0) = (*fluid)(i, j, 0);
		}
	}
	
	}
}

inline void calLatticeRhoVelocity()
{
    double rh=0.0, vx = 0.0, vy =0.0;


	for (size_t i = 1; i <= sizey; i++)
	{
		for (size_t j = 1; j <= sizex; j++)
		{
			rh = 0.0;
			vx = 0.0;
			vy = 0.0;
         if(flag && (*tmpfluid).getBoundary(i,j) == 0)
            {
			for (int k = 0; k < Q; k++)
			{
				rh += (*tmpfluid)(i, j, k);
				vx += (*tmpfluid)(i, j, k) * disvel[k][0];
				vy += (*tmpfluid)(i, j, k) * disvel[k][1];
			}

            (*tmpfluid)(i, j, 9) = vx ;
            (*tmpfluid)(i, j, 10) = vy;
			(*tmpfluid)(i, j, 11) = rh;
         }
		}
	}

}

inline double feq(size_t const k, double ux, double const uy, double const rho)
{
	double r = (ux * disvel[k][0] + uy * disvel[k][1]);
	double u2 = ux * ux + uy * uy;
	return (stencil[k] * (rho + 3.0 * r + 4.5 * r*r - 1.5 * u2));
}

inline void collide()
{
	for (size_t i = 1; i <= sizey; i++)
	{
		for (size_t j = 1; j <= sizex; j++)
		{
	 if(flag && (*tmpfluid).getBoundary(i,j) == 0)
     {
            for (int k = 0; k < Q; k++)
			{
			(*tmpfluid)(i, j, k) = (1.0 - omega) * (*tmpfluid)(i, j, k) + omega * feq(k, (*tmpfluid)(i, j, 9), (*tmpfluid)(i, j, 10), (*tmpfluid)(i, j, 11));
			}
     }
		}

	}

}


int main(int argc, char** argv)
{

	if (argc < 2)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	string fname = argv[1];
    string vtkfilename;
	size_t vtk_step;
	
	FileReader* fr = new FileReader();

	fr->readParameters( fname);
 	
	sizex = fr->getParameter<size_t>("sizex");
	sizey = fr->getParameter<size_t>("sizey");
	timesteps = fr->getParameter<size_t>("timesteps");
	omega = fr->getParameter<Real>("omega");
	vtkfilename = fr->getParameter<string>("vtk_file");
	vtk_step = fr->getParameter<size_t>("vtk_step");
	imfile = fr->getParameter<string>("geometry");

	std::cout << "sizex = " << sizex << '\n';
	std::cout << "sizey = " << sizey << '\n';
	std::cout << "timesteps = " << timesteps << '\n';
	std::cout << "omega = " << omega << '\n';
	std::cout << "vtk_file = " << vtkfilename << '\n';
	std::cout << "vtk_step = " << vtk_step << '\n';
	if(sizex == 9999999999 && imfile != "9999999999")
	flag = true;
	vtkfilename = vtkfilename.substr(0, vtkfilename.find('.'));
	std::cout << "vtk_file = " << vtkfilename << '\n';
    cout << flag << " ";
	init();

    timeval start, end;

    gettimeofday(&start, 0);



    for (size_t i = 0; i < timesteps; i++)
    {
        stream();
        calLatticeRhoVelocity();
        collide();
        (*fluid).copy(tmpfluid);
         if (i==1 || i%vtk_step == 0)
                {
                        string vtkfile = std::string("./output/" + vtkfilename) + std::string(to_string(i)) + std::string(".vtk");
                        writeVTK(vtkfile,fluid);
                       
                }
	}
    gettimeofday(&end, 0);
    double elapsed = 0.000001 * ((double)((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec));
    int nofluidcel = fluid->getFluidCellNumb();
    double mlu = (nofluidcel * timesteps) /elapsed ;
    std::cout << "\n MLUps Value = " << mlu << '\n';
	fluid->~Grid();
	tmpfluid->~Grid();
}
