#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace lbm {
	
	struct dat {
	  string key;
	  string value;
	};

	class FileReader
	{
		public:
			FileReader();
			void readParameters(string filename);
			
			template<typename Type>
			Type getParameter(string valueName)
			{
				for(vector<dat>::iterator it = m_params.begin() ; it != m_params.end(); ++it)
				{
					if(it->key == valueName)
					{
						return boost::lexical_cast<Type>(it->value);
					}
				}
				return boost::lexical_cast<Type>("nix wars");
			}
		
		private:
			vector<dat> m_params;
	};

}
