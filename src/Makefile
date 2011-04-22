LINKS=-L../lib

CXXFLAGS+= -O3 -g -Wall -fPIC -pg $(INCLUDES)
BINPATH=../bin

objects=output.o Dumper.o LatticeAnalysis.o \
	StatesClassification.o Term.o IndexClassification.o \
	HamiltonianPart.o Hamiltonian.o FieldOperatorPart.o FieldOperator.o \
	DensityMatrixPart.o DensityMatrix.o \
	GreensFunctionPart.o GreensFunction.o \
	TwoParticleGF.o TwoParticleGFPart.o \
	Vertex4.o FieldOperatorContainer.o \
	TwoParticleGFContainer.o GFContainer.o 

devel_objects= TwoParticleGFContainer.o

libs=../lib/libiniconfig.a ../lib/lib_json.a -liniparser # -lhdf5 -lhdf5_cpp

default: main
devel: $(devel_objects)

output.o: output.cpp output.h
Dumper.o: Dumper.cpp Dumper.h config.h
LatticeAnalysis.o: LatticeAnalysis.h LatticeAnalysis.cpp
Term.o:	Term.cpp Term.h
IndexClassification.o: IndexClassification.h IndexClassification.cpp Term.h
StatesClassification.o: StatesClassification.h StatesClassification.cpp config.h
HamiltonianPart.o: StatesClassification.h config.h HamiltonianPart.h HamiltonianPart.cpp
Hamiltonian.o: StatesClassification.h config.h HamiltonianPart.h Hamiltonian.h Hamiltonian.cpp
FieldOperatorPart.o: StatesClassification.h config.h HamiltonianPart.h FieldOperatorPart.h FieldOperatorPart.cpp
FieldOperator.o: StatesClassification.h config.h HamiltonianPart.h FieldOperatorPart.h FieldOperatorPart.cpp FieldOperator.h
FieldOperatorContainer.o : StatesClassification.o Hamiltonian.o FieldOperatorContainer.h FieldOperatorContainer.cpp 
DensityMatrixPart.o: DensityMatrixPart.h DensityMatrixPart.cpp
DensityMatrix.o: DensityMatrix.h DensityMatrix.cpp
GreensFunctionPart.o: GreensFunctionPart.h GreensFunctionPart.cpp config.h
GreensFunction.o: GreensFunction.h GreensFunction.cpp
TwoParticleGFPart.o: TwoParticleGFPart.h TwoParticleGFPart.cpp config.h
TwoParticleGF.o: TwoParticleGF.h TwoParticleGF.cpp
TwoParticleGFContainer.o : TwoParticleGFContainer.h TwoParticleGFContainer.cpp TwoParticleGFPart.o TwoParticleGF.o FieldOperatorContainer.o
GFContainer.o : GFContainer.h GFContainer.cpp GreensFunction.o GreensFunctionPart.o
Vertex4.o: Vertex4.h Vertex4.cpp TwoParticleGFContainer.o GFContainer.o

libpomerol.so: $(objects)
	$(CXX) -shared $(objects) $(libs) -o libpomerol.so

libpomerol.a: $(objects)
	ar -cvr libpomerol.a $(objects)

main : main.cpp libpomerol.a
	$(CXX) $(CXXFLAGS) main.cpp libpomerol.a $(libs) -o $(BINPATH)/pomerolDiag $(LINKS)

clean:
	$(RM) *.o *.a
