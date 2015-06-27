#include <cstdlib>
#include <string>
#include "Grid.h"
#include "typedefs.h"
#include "FileReader.h"
using namespace lbm;

int main(int args, char **argv)
{
	lbm::FileReader test;
	test.readParameters("referenceOutputs/params_5x5.dat");
	size_t timesteps = test.getParameter<size_t>("timesteps");
	size_t sizex = test.getParameter<size_t>("sizex");
	size_t sizey = test.getParameter<size_t>("sizey");
	
	printf("timesteps: %lu\n", timesteps);
	return 0;
}
