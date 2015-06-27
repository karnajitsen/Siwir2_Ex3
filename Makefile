CC=g++
CFLAGS= -Wall -std=c++11 -pedantic -g
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=lbm.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=lbm
COMMON=

all: clean lbm

lbm:
	$(CC) $(CFLAGS) $(SOURCES) -o lbm
	
prm:
	rm -rf ./output/*.*
	./lbm params.dat

ref0:
	./waveguide 0.01 0.0000000001 0
	rm -rf *.pdf
	gnuplot ./plot.p
test: 
	$(CC) $(CFLAGS) $(SOURCES) -o waveguide
	./waveguide 0.01 0.0000000001

clean:
	rm -rf lbm
	
.PHONY : all clean
