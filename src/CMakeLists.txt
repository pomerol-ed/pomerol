set(SOURCES
    mpi_dispatcher/mpi_dispatcher
    pomerol/Misc
    pomerol/LatticePresets
    pomerol/StatesClassification
    pomerol/HamiltonianPart
    pomerol/Hamiltonian
    pomerol/MonomialOperatorPart
    pomerol/MonomialOperator
    pomerol/FieldOperatorContainer
    pomerol/DensityMatrixPart
    pomerol/DensityMatrix
    pomerol/GreensFunctionPart
    pomerol/GreensFunction
    pomerol/GFContainer
    pomerol/TwoParticleGFPart
    pomerol/TwoParticleGF
    pomerol/TwoParticleGFContainer
    pomerol/Vertex4
    pomerol/SusceptibilityPart
    pomerol/Susceptibility
    pomerol/EnsembleAverage
)

add_library(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME}
                           PRIVATE
                           ${PROJECT_SOURCE_DIR}/include
                           ${PROJECT_BINARY_DIR}/include
                           PUBLIC
                           ${EIGEN3_INCLUDE_DIR}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${Boost_INCLUDE_DIRS})
set_property(TARGET ${PROJECT_NAME} APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC libcommute ${MPI_CXX_LIBRARIES})
target_compile_options(${PROJECT_NAME} PUBLIC ${MPI_CXX_COMPILE_FLAGS}
                                              ${MPI_C_COMPILE_FLAGS})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(OpenMP_CXX_FOUND)
    target_compile_options(${PROJECT_NAME} PRIVATE ${CMAKE_CXX_FLAGS}
                                                   ${OpenMP_CXX_FLAGS})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
endif()

# Install the main library

install(TARGETS ${PROJECT_NAME}
        EXPORT "${PROJECT_NAME}Targets"
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
