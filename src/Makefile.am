MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS= ${eigen3_CFLAGS} -I../jsoncpp/include/ -O3 -Wall -fPIC -g #-pg
lib_LIBRARIES=libpomerol.a
bin_PROGRAMS=pomerolDiag

libpomerol_a_SOURCES= output.cpp HDF5Storage.cpp Thermal.cpp LatticeAnalysis.cpp \
	StatesClassification.cpp Term.cpp IndexClassification.cpp \
	HamiltonianPart.cpp Hamiltonian.cpp FieldOperatorPart.cpp FieldOperator.cpp \
	DensityMatrixPart.cpp DensityMatrix.cpp \
	GreensFunctionPart.cpp GreensFunction.cpp \
	FourIndexObject.cpp \
	TwoParticleGF.cpp TwoParticleGFPart.cpp \
	Vertex4.cpp FieldOperatorContainer.cpp \
	TwoParticleGFContainer.cpp GFContainer.cpp Misc.h ComputableObject.h 
#libpomerol_a_HEADERS=Misc.h
pomerolDiag_SOURCES=main.cpp
pomerolDiag_LDFLAGS=${hdf5_LIBS} -lhdf5 -lhdf5_cpp libpomerol.a ../jsoncpp/src/libjsoncpp.a
pomerolDiag_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a

all-clean:
	$(RM) *.o *.a
