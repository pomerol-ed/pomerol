MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/

lib_LIBRARIES=libpomerol.a

libpomerol_a_CPPFLAGS=-Wall -g
libpomerol_a_SOURCES= \
	Lattice.h Lattice.cpp LatticePresets.h \
	Logger.cpp Logger.h 
#HDF5Storage.cpp Thermal.cpp Lattice.cpp \
#StatesClassification.cpp Term.cpp IndexClassification.cpp \
#HamiltonianPart.cpp Hamiltonian.cpp FieldOperatorPart.cpp FieldOperator.cpp \
#DensityMatrixPart.cpp DensityMatrix.cpp \
#GreensFunctionPart.cpp GreensFunction.cpp \
#FourIndexObject.cpp \
#TwoParticleGF.cpp TwoParticleGFPart.cpp \
#Vertex4.cpp FieldOperatorContainer.cpp \
#TwoParticleGFContainer.cpp GFContainer.cpp Misc.h ComputableObject.h Lattice.h
libpomerol_a_pkginclude_HEADERS=*.h 
libpomerol_a_pkgincludedir=$(includedir)/pomerol

all-clean:
	$(RM) *.o *.a
