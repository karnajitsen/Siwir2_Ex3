#include <iostream>
#include<fstream>
#include<cmath>
#include<memory>
#include<stdlib.h>

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

Grid *fluid, *tmpfluid;
double uw[2] = {0.8, 0};
double disvel[Q][2] = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 }, { -1.0, -1.0 }, { 1.0, -1.0 } };
int neighbours[Q][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };
double stencil[Q] = { 1.0 / 9.0, 4.0 / 9.0, 4.0 / 9.0, 4.0 / 9.0, 4.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

inline void init()
{
cout << "@@@@################";
	double ux = 0.0, uy = 0.0, initrho =1.0;
	fluid = new Grid(sizex, sizey, stencil, ux, uy, initrho);
	tmpfluid = new Grid(sizex, sizey, stencil, ux, uy, initrho);
	/*for (int i = 0; i < Q; i++)
	feq[i] = stencil[i];*/
cout << "@@@@@@@@@@@@@@@@@@@@@@";
}

inline double calLidVel(size_t k)
{
	return (6.0 * stencil[k] * (disvel[k][0] * uw[0] + disvel[k][1] * uw[1]));
}

inline void stream()
{
	for (size_t i = 1; i <= sizey; i++)
	{
		for (size_t j = 1; j <= sizex; j++)
		{
			for (int k = 1; k < Q; k++)
			{
				double pt = (*tmpfluid)(i + neighbours[k][0], j + neighbours[k][1], 12);
				size_t l;
                                        if(k != 4)
                                        l = (k+((Q-1)/2))% ( Q-1);
                                        else
                                        l=8;

				if (pt > 0.0 && pt != 2.0 )
					(*tmpfluid)(i, j, l) = (*fluid)(i, j, k);       // Bounce back from vertical and bottom walls
					
				if (pt == 2.0)
					(*tmpfluid)(i, j, l) = (*fluid)(i, j, k) - calLidVel(k); // Streaming effect due to lid velocity in upper wall
				else
				    (*tmpfluid)(i + neighbours[k][0], j + neighbours[k][1], k) = (*fluid)(i,j,k);  // Free Streaming
				
			}

			(*tmpfluid)(i, j, 0) = (*fluid)(i, j, 0);
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
			for (int k = 0; k < Q; k++)
			{
				rh += (*tmpfluid)(i, j, k);
				vx += (*tmpfluid)(i, j, k) * disvel[k][0];
				vy += (*tmpfluid)(i, j, k) * disvel[k][1];
			}

			(*tmpfluid)(i, j, 9) = vx / rh;
			(*tmpfluid)(i, j, 10) = vy / rh;
			(*tmpfluid)(i, j, 11) = rh;

		}
	}

}

inline double feq(size_t const k, double ux, double const uy, double const rho)
{
	double r = (ux * disvel[k][0] + uy * disvel[k][1]);
	double u2 = ux * ux + uy * uy;
	return stencil[k] * (rho + 3.0 * r + 4.5 * r*r - 1.5 * u2);
}

inline void collide()
{
	for (size_t i = 1; i <= sizey; i++)
	{
		for (size_t j = 1; j <= sizex; j++)
		{
			for (int k = 0; k < Q; k++)
			{
				(*tmpfluid)(i, j, k) = (1.0 - omega) * (*tmpfluid)(i, j, k) + feq(k, (*tmpfluid)(i, j, 9), (*tmpfluid)(i, j, 10), (*tmpfluid)(i, j, 11));
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
	ifstream paramfile;
	string tmp;
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

	std::cout << "sizex = " << sizex << '\n';
	std::cout << "sizey = " << sizey << '\n';
	std::cout << "timesteps = " << timesteps << '\n';
	std::cout << "omega = " << omega << '\n';
	std::cout << "vtk_file = " << vtkfilename << '\n';
	std::cout << "vtk_step = " << vtk_step << '\n';

	vtkfilename = vtkfilename.substr(0, vtkfilename.find('.'));
	std::cout << "vtk_file = " << vtkfilename << '\n';

	cout << "$$$$$$$$$$";
	init();
	int k = 0;
	
	for (size_t i = 0; i < timesteps; i++)
	{
		if (i%vtk_step == 0)
		{
			string vtkfile = std::string("./output/" + vtkfilename) + std::string(to_string(k)) + std::string(".vtk");
			writeVTK(vtkfile,(*fluid));
		}
		stream();
		calLatticeRhoVelocity();
		collide();
		(*fluid).copy((*tmpfluid));	

	}
	fluid->~Grid();
	tmpfluid->~Grid();
}
