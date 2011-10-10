SOURCES = src/main.cpp ../amplitudelib/pdf/mrst99.cpp ../amplitudelib/tools/tools.cpp ../amplitudelib/amplitudelib/amplitudelib.cpp \
		../amplitudelib/amplitudelib/datafile.cpp ../amplitudelib/tools/interpolation.cpp \
		src/fourpoint_fft.cpp ../amplitudelib/pdf/pdf.cpp \
		../amplitudelib/pdf/mrst.cpp ../amplitudelib/pdf/cteq.cpp src/xs.cpp \
		../amplitudelib/fragmentation/fragmentation.cpp ../amplitudelib/fragmentation/kkp.cpp \
		src/fourpoint.cpp ../amplitudelib/amplitudelib/virtual_photon.cpp \
		../amplitudelib/amplitudelib/wave_function.cpp
CSOURCES = ../amplitudelib/fourier/fourier.c
FSOURCES = ../amplitudelib/fragmentation/fragmentation_kkp.f ../amplitudelib/pdf/CT10Pdf.f
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.c=.o)
FOBJECTS=$(FSOURCES:.f=.o)


