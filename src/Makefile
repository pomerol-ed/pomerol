LINKS=-L../lib

CXXFLAGS+= -O3 -g -Wall -fPIC -pg $(INCLUDES)
BINPATH=../bin

objects=output.o LatticeAnalysis.o \
	StatesClassification.o Term.o BitClassification.o \
	HamiltonianPart.o Hamiltonian.o FieldOperatorPart.o FieldOperator.o \
	DensityMatrixPart.o DensityMatrix.o \
	GreensFunctionPart.o GreensFunction.o \
	TwoParticleGF.o TwoParticleGFPart.o \
	Vertex4.o ContainerTwoParticleGF.o

devel_objects= ContainerTwoParticleGF.o

libs=../lib/libiniconfig.a ../lib/lib_json.a -liniparser

default: main
devel: $(devel_objects)

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
ContainerTwoParticleGF.o : ContainerTwoParticleGF.h ContainerTwoParticleGF.cpp
Vertex4.o: Vertex4.h Vertex4.cpp

libmeat.so: $(objects)
	$(CXX) -shared $(objects) $(libs) -o libmeat.so

libmeat.a: $(objects)
	ar -cvr libmeat.a $(objects)

main : main.cpp libmeat.a
	$(CXX) $(CXXFLAGS) main.cpp libmeat.a $(libs) -o $(BINPATH)/pomerolDiag $(LINKS)

clean:
	$(RM) *.o *.a
