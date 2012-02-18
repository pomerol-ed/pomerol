MAINTAINERCLEANFILES = Makefile.in
AM_CXXFLAGS+= ${hdf5_CFLAGS} ${eigen3_CFLAGS} -I../jsoncpp/include/

lib_LIBRARIES=libpomerol.a

libpomerol_a_CPPFLAGS=-Wall -g
libpomerol_a_SOURCES= Misc.cpp Misc.h \
	HDF5Storage.cpp HDF5Storage.h \
	Logger.cpp Logger.h \
	Thermal.cpp Thermal.h ComputableObject.h \
	Index.cpp Index.h \
	LatticeReader.cpp LatticeReader.h \
	LatticeAnalysis.cpp LatticeAnalysis.h \
	Term.cpp Term.h \
	IndexClassification.cpp IndexClassification.h \
	StatesClassification.cpp StatesClassification.h \
	HamiltonianPart.cpp HamiltonianPart.h Hamiltonian.cpp Hamiltonian.h \
	FieldOperatorPart.cpp FieldOperatorPart.h FieldOperator.cpp FieldOperator.h \
	FieldOperatorContainer.cpp FieldOperatorContainer.h \
	DensityMatrixPart.cpp DensityMatrixPart.h DensityMatrix.cpp DensityMatrix.h \
	IndexContainer2.h \
	GreensFunctionPart.cpp GreensFunctionPart.h GreensFunction.cpp GreensFunction.h \
	GFContainer.cpp GFContainer.h \
	IndexContainer4.h \
	TwoParticleGFPart.cpp TwoParticleGFPart.h TwoParticleGF.cpp TwoParticleGF.h \
	TwoParticleGFContainer.cpp TwoParticleGFContainer.h
#	Vertex4.cpp Vertex4.h \
#	Vertex4Container.cpp Vertex4Container.h \
#	AmputatedVertex4.cpp AmputatedVertex4.h \
#	AmputatedVertex4Container.cpp AmputatedVertex4Container.h
	
libpomerol_a_pkginclude_HEADERS=*.h
libpomerol_a_pkgincludedir=$(includedir)/pomerol

all-clean:
	$(RM) *.o *.a
