include_directories(${CMAKE_SOURCE_DIR}/include/pomerol)

set (pomerol_src Misc.cpp 
	Lattice.cpp 
    LatticePresets.cpp 
	Index.cpp 
	IndexClassification.cpp 
	Operator.cpp 
    OperatorPresets.cpp 
	IndexHamiltonian.cpp 
	Symmetrizer.cpp 
	StatesClassification.cpp 
	#HDF5Storage.h HDF5Storage.cpp 
	HamiltonianPart.cpp 
 	Hamiltonian.cpp 
	Thermal.cpp 
	FieldOperatorPart.cpp 
    FieldOperator.cpp 
    FieldOperatorContainer.cpp 
	DensityMatrixPart.cpp 
    DensityMatrix.cpp 
	GreensFunctionPart.cpp 
    GreensFunction.cpp
	GFContainer.cpp
	TwoParticleGFPart.cpp 
    TwoParticleGF.cpp
	TwoParticleGFContainer.cpp 
	Vertex4.cpp 
	Vertex4Container.cpp 
    $<TARGET_OBJECTS:mpi_dispatcher>
)
add_library(pomerol SHARED ${pomerol_src})
set_target_properties(pomerol PROPERTIES VERSION ${POMEROL_VERSION})

target_link_libraries(pomerol
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    #${HDF5_LIBRARIES}
)
set_target_properties(pomerol PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS pomerol 
        EXPORT pomerol
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
)
install(EXPORT pomerol NAMESPACE pomerol:: DESTINATION share/pomerol)

