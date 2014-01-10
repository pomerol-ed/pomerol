include_directories (.)

set (pomerol_src Misc.cpp 
	Lattice.h Lattice.cpp LatticePresets.h LatticePresets.cpp 
    MPIDispatcher.h MPIDispatcher.cpp
	Index.h Index.cpp 
	IndexClassification.h IndexClassification.cpp 
	Operator.h Operator.cpp OperatorPresets.h OperatorPresets.cpp 
	IndexHamiltonian.h IndexHamiltonian.cpp 
	Symmetrizer.h Symmetrizer.cpp 
	StatesClassification.h StatesClassification.cpp 
	#HDF5Storage.h HDF5Storage.cpp 
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
set_target_properties(pomerol PROPERTIES VERSION ${POMEROL_VERSION})

target_link_libraries(pomerol
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    #${HDF5_LIBRARIES}
)

install(TARGETS pomerol
    DESTINATION lib
    )

install(DIRECTORY . DESTINATION include/pomerol
          FILES_MATCHING PATTERN "*.h")
