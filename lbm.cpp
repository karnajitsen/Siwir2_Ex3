#include <iostream>
#include<fstream>
#include<cmath>
#include<memory>
#include<stdlib.h>

#include "Grid.h"
typedef double Real;


//Global Variables
size_t sizex, sizey, timesteps;
Real omega;


int main(int argc, char** argv)
{

	if (argc < 4)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	string fname = argv[1];
	ifstream paramfile;
	string tmp;
	paramfile.open("./LBMreferenceOutputs/referenceOutputs/" + fname);
	
	paramfile >> tmp;
	paramfile >> sizex;

	paramfile >> tmp;
	paramfile >> sizey;

	paramfile >> tmp;
	paramfile >> timesteps;

	paramfile >> tmp;
	paramfile >> omega;



}