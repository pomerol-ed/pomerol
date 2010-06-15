JSONCPP_PATH=../jsoncpp/include
EIGEN_INCLUDE=$(shell pkg-config --cflags-only-I eigen2)
INCLUDES=-I$(JSONCPP_PATH) $(EIGEN_INCLUDE)

LINKS=-L../lib
# FIXME
ifeq ($(shell uname),Darwin)
LINKS += -L/opt/local/lib ../lib/libiniparser.a
INCLUDES += -I/opt/local/include
else
LINKS += -liniparser
ifneq ($(shell echo $(EIGEN_INCLUDE) | head -c2),"-I")
	INCLUDES += -I/usr/include/eigen2
	EIGEN_INCLUDE=
endif
endif

CXXFLAGS= -O3 -g -Wall -fPIC -pg $(INCLUDES)
BINPATH=../bin

objects=iniconfig.o output.o LatticeAnalysis.o \
	StatesClassification.o Term.o BitClassification.o \
	HamiltonianPart.o Hamiltonian.o FieldOperatorPart.o FieldOperator.o \
	DensityMatrixPart.o DensityMatrix.o \
	GreensFunctionPart.o GreensFunction.o \
	TwoParticleGF.o TwoParticleGFPart.o

default :  main

iniconfig.o: iniconfig.cpp iniconfig.h
output.o: output.cpp output.h
LatticeAnalysis.o: LatticeAnalysis.h LatticeAnalysis.cpp
Term.o:	Term.cpp Term.h
BitClassification.o: BitClassification.h BitClassification.cpp Term.h
StatesClassification.o: StatesClassification.h StatesClassification.cpp config.h
HamiltonianPart.o: StatesClassification.h config.h HamiltonianPart.h HamiltonianPart.cpp
Hamiltonian.o: StatesClassification.h config.h HamiltonianPart.h Hamiltonian.h Hamiltonian.cpp
FieldOperatorPart.o: StatesClassification.h config.h HamiltonianPart.h FieldOperatorPart.h FieldOperatorPart.cpp
FieldOperator.o: StatesClassification.h config.h HamiltonianPart.h FieldOperatorPart.h FieldOperatorPart.cpp FieldOperator.h
DensityMatrixPart.o: DensityMatrixPart.h DensityMatrixPart.cpp
DensityMatrix.o: DensityMatrix.h DensityMatrix.cpp
GreensFunctionPart.o: GreensFunctionPart.h GreensFunctionPart.cpp config.h
GreensFunction.o: GreensFunction.h GreensFunction.cpp
TwoParticleGFPart.o: TwoParticleGFPart.h TwoParticleGFPart.cpp config.h
TwoParticleGF.o: TwoParticleGF.h TwoParticleGF.cpp 

libmeat.so: $(objects)
	$(CXX) -shared $(objects) -o libmeat.so

libmeat.a: $(objects)
	ar -cvr libmeat.a $(objects)

main : main.cpp libmeat.a
	$(CXX) $(CXXFLAGS) main.cpp libmeat.a ../lib/lib_json.a -o $(BINPATH)/pomerolDiag $(LINKS)

clean:
	rm -f *.o *.a
