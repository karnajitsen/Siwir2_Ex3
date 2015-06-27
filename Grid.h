#pragma once
#include<iostream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string.h>
#include<malloc.h>
#define LD 16
#define ALLIGNMENT 32
#define CELLS 14
using namespace std;
//#define M_PI 3.14
class Grid
{

    //__declspec(align(128))
    double * __restrict data = NULL;
    size_t sizeX, sizeY, ld, size;
    
public:
    explicit Grid()
    {
        data = (double*)memalign(ALLIGNMENT, 0);
        sizeX = 0;
        sizeY = 0;
		size = 0.0;
    }

	explicit Grid(const size_t x, const size_t y, const double* fi, const double ux, const double uy, const double rho)
    {
        sizeX = x+2;
        sizeY = y+2;
		ld = sizeX + LD;
		size = CELLS*ld*y*sizeof(double);
        //data = (double*) memalign(ALLIGNMENT, CELLS*ld*y*sizeof(double));
		data = new double[CELLS*ld*y];
        //data = (double*) _aligned_malloc(ld*y*sizeof(double), ALLIGNMENT);
		for (size_t i = 0; i < sizeY; i++)
		{
			for (size_t j = 0; j < sizeY; j++)
			{
				
				data[i * ld + j * CELLS + 0] = fi[0];
				data[i * ld + j * CELLS + 1] = fi[1];
				data[i * ld + j * CELLS + 2] = fi[2];
				data[i * ld + j * CELLS + 3] = fi[3];
				data[i * ld + j * CELLS + 4] = fi[4];
				data[i * ld + j * CELLS + 5] = fi[5];
				data[i * ld + j * CELLS + 6] = fi[6];
				data[i * ld + j * CELLS + 7] = fi[7];
				data[i * ld + j * CELLS + 8] = fi[8];
				data[i * ld + j * CELLS + 9] = ux;
				data[i * ld + j * CELLS + 10] = uy;
				data[i * ld + j * CELLS+12] = 0.0;
				if (i == 0) data[i * ld + j * CELLS+12] = 1.0;
				if (j == 0) data[i * ld + j * CELLS+12] = 4.0;
				if (j == sizeX - 1)  data[i * ld + j * CELLS+12] = 2.0;
				if (i == sizeY - 1) data[i * ld + j * CELLS+12] = 3.0;
				data[i * ld + j * CELLS + 11] =  rho;
			}
		}

        //data++;
    }
    ~Grid()
    {
        //--data;
        free(data);
    }

	inline void copy(Grid grd)
	{
		memcpy((void *)data, (void *)&grd(0, 0, 0), size);
	}

    /*inline void reset()
    {
        for (size_t i = 0; i < sizeY; i++)
        {
            for (size_t j = 0; j < sizeX; j++)
            {
                data[i*ld + j] = 0.0;
            }
        }
    }*/

	inline void swap(const Grid tgtgrd)
	{

	}

    inline double& operator()(const size_t x, const size_t y, const size_t f)
    {
        assert(x < sizeX);
        assert(y < sizeY);
        return data[y*ld + x * CELLS + f];
    }

	inline double& operator()(const size_t x, const size_t y, const size_t f) const
    {
        assert(x < sizeX);
        assert(y < sizeY);
		return data[y*ld + x * CELLS + f];
    }

    inline Grid * operator+=(const Grid * rhs)
    {
        return this;
    }

    inline Grid* operator-=(const Grid* rhs)
    {
        return this;
    }

    inline size_t getXsize() const
    {
        return sizeX;
    }

    inline size_t getYsize() const
    {
        return sizeY;
    }

   
};
