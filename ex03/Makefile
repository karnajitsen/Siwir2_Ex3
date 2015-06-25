TARGET = lbm
CC     = g++
#FLAGS  = -O3 -Wall -Winline -Wshadow -std=c++11
FLAGS  = -pg -Wall -Wshadow -std=c++11
.PHONY: clean
default: clean lbm

lbm: main.cpp
	$(CC) $(FLAGS) main.cpp Grid.cpp FileReader.cpp -o $(TARGET)
	
clean:
	rm $(TARGET)
	

