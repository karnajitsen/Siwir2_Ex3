CC=g++
CFLAGS= -Wall -std=c++11 -pedantic -g
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=waveguide.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=waveguide
COMMON=

all: wave

wave:
	$(CC) $(CFLAGS) $(SOURCES) -o waveguide
	
ref1:
	./waveguide 0.01 0.0000000001 1

ref0:
	./waveguide 0.01 0.0000000001 0
	rm -rf *.pdf
	gnuplot ./plot.p
test: 
	$(CC) $(CFLAGS) $(SOURCES) -o waveguide
	./waveguide 0.01 0.0000000001

clean:
	
.PHONY : all clean
