SOURCES = src/main.cpp src/pdf/mrst99.cpp ../amplitudelib/tools/tools.cpp ../amplitudelib/amplitudelib/amplitudelib.cpp \
		../amplitudelib/amplitudelib/datafile.cpp ../amplitudelib/tools/interpolation.cpp \
		 src/pdf.cpp \
		src/pdf/mrst.cpp src/pdf/cteq.cpp src/xs.cpp \
		src/fragmentation/fragmentation.cpp src/fragmentation/kkp.cpp \
		src/fourpoint.cpp
CSOURCES = ../amplitudelib/fourier/fourier.c
FSOURCES = src/fragmentation/fragmentation_kkp.f src/pdf/CT10Pdf.f
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.c=.o)
FOBJECTS=$(FSOURCES:.f=.o)


