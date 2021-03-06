#pragma once

#include <cstdlib>
#include <utility>

namespace lbm {

	template<typename Type, size_t Cellsize>
	class Grid
	{
		public:
			inline Grid();
			inline Grid(size_t xsize, size_t ysize);
			inline ~Grid();
			inline Type& operator()(size_t x, size_t y, size_t f);
			inline Type  operator()(size_t x, size_t y, size_t f) const;
			inline void swap( Grid& grid );
			
		private:
			size_t m_sizeX;
			size_t m_sizeY;
			Type* m_data;
	};
	
	
	template< typename Type >
	class Grid<Type,1>
	{
		public:
			inline Type& operator()( size_t x, size_t y );
			inline Type  operator()( size_t x, size_t y ) const;
	};
	
	
	template< typename Type >
	class Grid<Type,0>;
	
}
