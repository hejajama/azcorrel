CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp -I ../amplitudelib# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm -lfftw3 -lgfortran -lpthread #-L../amplitudelib/ -lamplitude 
CFLAGS = -O2 -D_REENTRANT

#CC = gcc
#CXX = g++
#FOR = gfortran
CC = mpicc
CXX = mpic++
FOR = mpif90

include filelist.m

all: azcorrel 

azcorrel: $(OBJECTS) $(COBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(COBJECTS) ../amplitudelib/libamplitude.a -lgfortran -o azcorrel
.cpp.o:
	 $(CXX) $(CXXFLAGS) $< -c -o $@
.c.o:
	$(CC) $(CXXFLAGS) $(CFLAGS) $< -c -o $@
.f.o:
	$(FOR) -c $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(COBJECTS)
	rm -f azcorrel
