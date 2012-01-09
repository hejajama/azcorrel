#SOURCES = src/main.cpp ../amplitudelib/pdf/mrst99.cpp ../amplitudelib/tools/tools.cpp ../amplitudelib/amplitudelib/amplitudelib.cpp \
#		../amplitudelib/amplitudelib/datafile.cpp ../amplitudelib/tools/interpolation.cpp \
#		src/fourpoint_fft.cpp ../amplitudelib/pdf/pdf.cpp \
#		../amplitudelib/pdf/mrst.cpp ../amplitudelib/pdf/cteq.cpp src/xs.cpp \
#		../amplitudelib/fragmentation/fragmentation.cpp ../amplitudelib/fragmentation/kkp.cpp \
#		src/fourpoint.cpp ../amplitudelib/amplitudelib/virtual_photon.cpp \
#		../amplitudelib/amplitudelib/wave_function.cpp ../amplitudelib/tools/interpolation2d.cpp
SOURCES = src/main.cpp src/fourpoint_fft.cpp src/fourpoint.cpp src/xs.cpp 
CSOURCES = src/pvegas/pvegas_mpi.c
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.c=.o)

