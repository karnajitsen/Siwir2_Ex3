#include "Grid.h"

namespace lbm{

	template<typename Type, size_t Cellsize>
	Grid<Type,Cellsize>::Grid()
	: m_sizeX(0), m_sizeY(0), m_data(0)
	{
	
	}

	template<typename Type, size_t Cellsize>
	Grid<Type,Cellsize>::Grid(size_t xsize, size_t ysize)
	: m_sizeX(xsize), m_sizeY(ysize), m_data( new Type[Cellsize*xsize*ysize] )
	{
	
	}


	template<typename Type, size_t Cellsize>
	inline Type& Grid<Type,Cellsize>::operator()(size_t x, size_t y, size_t f)
	{
		this->assert(x < m_sizeX && y < m_sizeY && f < Cellsize);
		return m_data[y*m_sizeX*Cellsize+x*Cellsize+f];
	}
	

	template<typename Type, size_t Cellsize>
	inline Type Grid<Type,Cellsize>::operator()(size_t x, size_t y, size_t f) const
	{
		this->assert(x < m_sizeX && y < m_sizeY && f < Cellsize);
		return m_data[y*m_sizeX*Cellsize+x*Cellsize+f];
	}


	template<typename Type, size_t Cellsize>
	inline void Grid<Type,Cellsize>::swap(Grid& grid)
	{
		std::swap(m_sizeY, grid.m_sizeY);
		std::swap(m_sizeX, grid.m_sizeX);
		std::swap(m_data, grid.m_data);
	}
	
	
	// Global swap function
	template<typename Type, size_t N>
	inline void swap(Grid<Type,N>& a, Grid<Type,N>& b)
	{
		a.swap(b);
	}

}
