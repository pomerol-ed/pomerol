LINKS=-L../lib

CXXFLAGS+= -O3 -Wall -fPIC -g $(INCLUDES) #-pg
BINPATH=../bin

objects=output.o Dumper.o LatticeAnalysis.o \
	StatesClassification.o Term.o IndexClassification.o \
	HamiltonianPart.o Hamiltonian.o FieldOperatorPart.o FieldOperator.o \
	DensityMatrixPart.o DensityMatrix.o \
	GreensFunctionPart.o GreensFunction.o \
	FourIndexObject.o \
	TwoParticleGF.o TwoParticleGFPart.o \
	Vertex4.o FieldOperatorContainer.o \
	TwoParticleGFContainer.o GFContainer.o 

libs=../lib/lib_json.a ${POMEROL_LIBS}

default: main

output.o: output.cpp output.h
Dumper.o: Dumper.cpp Dumper.h Misc.h
LatticeAnalysis.o: LatticeAnalysis.h LatticeAnalysis.cpp
Term.o:	Term.cpp Term.h
IndexClassification.o: IndexClassification.h IndexClassification.cpp Term.h
StatesClassification.o: StatesClassification.h StatesClassification.cpp Misc.h
HamiltonianPart.o: StatesClassification.h Misc.h HamiltonianPart.h HamiltonianPart.cpp ComputableObject.h
Hamiltonian.o: StatesClassification.h Misc.h HamiltonianPart.h Hamiltonian.h Hamiltonian.cpp ComputableObject.h
FieldOperatorPart.o: StatesClassification.h Misc.h HamiltonianPart.h FieldOperatorPart.h FieldOperatorPart.cpp ComputableObject.h
FieldOperator.o: StatesClassification.h Misc.h HamiltonianPart.h FieldOperatorPart.h FieldOperator.h FieldOperator.cpp ComputableObject.h
FieldOperatorContainer.o : StatesClassification.h Hamiltonian.h FieldOperatorContainer.h FieldOperatorContainer.cpp ComputableObject.h
DensityMatrixPart.o: DensityMatrixPart.h DensityMatrixPart.cpp ComputableObject.h
DensityMatrix.o: DensityMatrix.h DensityMatrix.cpp ComputableObject.h
GreensFunctionPart.o: GreensFunctionPart.h GreensFunctionPart.cpp Misc.h ComputableObject.h
GreensFunction.o: GreensFunction.h GreensFunction.cpp ComputableObject.h
TwoParticleGFPart.o: TwoParticleGFPart.h TwoParticleGFPart.cpp Misc.h ComputableObject.h
FourIndexObject.o: FourIndexObject.h FourIndexObject.cpp Misc.h
TwoParticleGF.o: TwoParticleGF.h TwoParticleGF.cpp FourIndexObject.h ComputableObject.h
TwoParticleGFContainer.o : TwoParticleGFContainer.h TwoParticleGFContainer.cpp FourIndexObject.h TwoParticleGFPart.h ComputableObject.h
GFContainer.o : GFContainer.h GFContainer.cpp GreensFunction.h GreensFunctionPart.h ComputableObject.h
Vertex4.o: Vertex4.h Vertex4.cpp TwoParticleGFContainer.h GFContainer.h FourIndexObject.h ComputableObject.h

libpomerol.so: $(objects)
	$(CXX) -shared $(objects) $(libs) -o libpomerol.so

libpomerol.a: $(objects)
	ar -cvr libpomerol.a $(objects)

main : main.cpp libpomerol.a
	$(CXX) $(CXXFLAGS) main.cpp libpomerol.a $(libs) -o $(BINPATH)/pomerolDiag $(LINKS)

clean:
	$(RM) *.o *.a
