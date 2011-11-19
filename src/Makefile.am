MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/

lib_LIBRARIES=libpomerol.a

libpomerol_a_CPPFLAGS=-Wall -g
libpomerol_a_SOURCES= HDF5Storage.cpp Thermal.cpp LatticeAnalysis.cpp \
	StatesClassification.cpp Term.cpp IndexClassification.cpp \
	HamiltonianPart.cpp Hamiltonian.cpp FieldOperatorPart.cpp FieldOperator.cpp \
	DensityMatrixPart.cpp DensityMatrix.cpp \
	GreensFunctionPart.cpp GreensFunction.cpp \
	FourIndexObject.cpp \
	TwoParticleGF.cpp TwoParticleGFPart.cpp \
	Vertex4.cpp FieldOperatorContainer.cpp \
	TwoParticleGFContainer.cpp GFContainer.cpp Misc.h ComputableObject.h LatticeAnalysis.h
libpomerol_a_include_HEADERS=*.h 
libpomerol_a_includedir=$(includedir)/pomerol

all-clean:
	$(RM) *.o *.a
