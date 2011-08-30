CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp -I ../amplitudelib# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm -lfftw3 

include filelist.m

all: azcorrel 

azcorrel: $(OBJECTS) $(COBJECTS) $(FOBJECTS) 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(COBJECTS) $(FOBJECTS) -lgfortran -o azcorrel
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
.f.o:
	gfortran -c $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(COBJECTS)
	rm -f $(FOBJECTS)
	rm -f azcorrel
