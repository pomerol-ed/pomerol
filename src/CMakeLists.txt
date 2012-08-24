include_directories (.)

set (pomerol_src Misc.cpp 
	Lattice.h Lattice.cpp LatticePresets.h LatticePresets.cpp 
	Logger.cpp Logger.h 
	Index.h Index.cpp 
	IndexClassification.h IndexClassification.cpp 
	Operator.h Operator.cpp OperatorPresets.h OperatorPresets.cpp 
	IndexHamiltonian.h IndexHamiltonian.cpp 
	Symmetrizer.h Symmetrizer.cpp 
	StatesClassification.h StatesClassification.cpp 
	HDF5Storage.h HDF5Storage.cpp 
	HamiltonianPart.cpp 
 	Hamiltonian.cpp 
	Thermal.cpp Thermal.h 
	FieldOperatorPart.cpp FieldOperator.cpp FieldOperatorContainer.cpp 
	DensityMatrixPart.cpp DensityMatrixPart.h DensityMatrix.cpp DensityMatrix.h 
	IndexContainer2.h 
	GreensFunctionPart.cpp GreensFunctionPart.h GreensFunction.cpp GreensFunction.h 
	GFContainer.cpp GFContainer.h 
	IndexContainer4.h 
	TwoParticleGFPart.cpp TwoParticleGFPart.h TwoParticleGF.cpp TwoParticleGF.h 
	TwoParticleGFContainer.cpp TwoParticleGFContainer.h 
	Vertex4.cpp Vertex4.h 
	Vertex4Container.cpp Vertex4Container.h
)
add_library(pomerol SHARED ${pomerol_src})

target_link_libraries(pomerol
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    jsoncpp
)

install(TARGETS pomerol
    LIBRARY DESTINATION lib
    )

install(DIRECTORY . DESTINATION include/pomerol
          FILES_MATCHING PATTERN "*.h")
