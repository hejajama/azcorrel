CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp -I ../amplitudelib# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm -lfftw3 -lgfortran -lpthread #-L../amplitudelib/ -lamplitude 
CFLAGS = -O2 -D_REENTRANT

include filelist.m

all: azcorrel 

azcorrel: $(OBJECTS) $(COBJECTS) 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(COBJECTS) ../amplitudelib/libamplitude.a -lgfortran -o azcorrel
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $(CFLAGS) $< -c -o $@
.f.o:
	gfortran -c $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(COBJECTS)
	rm -f azcorrel
