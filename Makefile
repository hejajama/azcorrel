CXXFLAGS = `gsl-config --cflags` -Wno-long-long -O2 -Wall -pedantic -I ../amplitudelib #-fopenmp 
LDFLAGS = `gsl-config --libs` -lm -lpthread -lgfortran #-L../amplitudelib/ -lamplitude 
CFLAGS = -D_REENTRANT -std=c99
FORTRANFLAGS = 

#CC = gcc
#CXX = g++
#FOR = gfortran
CC = mpicc
CXX = mpic++
FOR = mpif90

include filelist.m

all: azcorrel 

azcorrel: $(OBJECTS) $(COBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(COBJECTS) ../amplitudelib/libamplitude_mpi.a -o azcorrel
.cpp.o:
	 $(CXX) $(CXXFLAGS) $< -c -o $@
.c.o:
	$(CC) $(CXXFLAGS) $(CFLAGS) $< -c -o $@
.f.o:
	$(FOR) $(FORTRANFLAGS) -c $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(COBJECTS)
	rm -f azcorrel
