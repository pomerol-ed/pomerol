MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS= ${AMM_CXXFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/ #${pomerol_CFLAGS} -O3 #-Wall -fPIC -g $(INCLUDES) #-pg
lib_LIBRARIES=libpomerol.a
bin_PROGRAMS=pomerolDiag

libpomerol_a_SOURCES= output.cpp Dumper.cpp LatticeAnalysis.cpp \
	StatesClassification.cpp Term.cpp IndexClassification.cpp \
	HamiltonianPart.cpp Hamiltonian.cpp FieldOperatorPart.cpp FieldOperator.cpp \
	DensityMatrixPart.cpp DensityMatrix.cpp \
	GreensFunctionPart.cpp GreensFunction.cpp \
	FourIndexObject.cpp \
	TwoParticleGF.cpp TwoParticleGFPart.cpp \
	Vertex4.cpp FieldOperatorContainer.cpp \
	TwoParticleGFContainer.cpp GFContainer.cpp Misc.h ComputableObject.h 
pomerolDiag_SOURCES=main.cpp
pomerolDiag_LDFLAGS=${hdf5_LIBS} libpomerol.a ../jsoncpp/src/libjsoncpp.a
pomerolDiag_DEPENDENCIES = $(top_builddir)/jsoncpp/src/libjsoncpp.a

all-clean:
	$(RM) *.o *.a
