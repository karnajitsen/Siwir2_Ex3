#include "FileReader.h"

namespace lbm {
	
	FileReader::FileReader()
	{
		
	}
	
	void FileReader::readParameters(string filename)
	{
		ifstream input;
		string line;
	
		input.open(filename, ios::in);
		if(input.is_open())
		{
			while(getline(input,line))
			{
				size_t space = line.find_first_of(' ');
				dat param;
				param.key = line.substr(0, space);
				param.value = line.substr(space+1);
				m_params.push_back(param);
			}
			input.close();
		}
	}

}
