SOURCES = src/main.cpp src/pdf/mrst99.cpp src/tools.cpp src/amplitudelib/amplitudelib.cpp \
		src/amplitudelib/datafile.cpp src/interpolation.cpp src/pdf.cpp \
		src/pdf/mrst.cpp src/pdf/cteq.cpp src/xs.cpp \
		src/fragmentation/fragmentation.cpp src/fragmentation/kkp.cpp
CSOURCES = src/amplitudelib/fourier/fourier.c
FSOURCES = src/fragmentation/fragmentation_kkp.f src/pdf/CT10Pdf.f
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.c=.o)
FOBJECTS=$(FSOURCES:.f=.o)


