cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project (pomerol)
set(POMEROL_VERSION 2.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
 set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Select default build type
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

#Use complex Matrix elements
option(use_complex "Use complex matrix elements" OFF)
if (use_complex)
    message(STATUS "Using complex matrix elements")
    add_definitions("-DPOMEROL_COMPLEX_MATRIX_ELEMENS")
endif (use_complex)

# Documentation
set(DOXYFILE_SOURCE_DIR "./src")
include(UseDoxygen)

#set(Boost_USE_STATIC_LIBS        ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)
find_package (MPI)
    message(STATUS "MPI includes: " ${MPI_INCLUDE_DIRS} )
    message(STATUS "MPI CXX libs: " ${MPI_CXX_LIBRARIES} )
    message(STATUS "MPI definitions : " ${MPI_CXX_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
    add_definitions(${MPI_CXX_COMPILE_FLAGS}  ${MPI_C_COMPILE_FLAGS})

find_package (Boost COMPONENTS mpi serialization)
    message(STATUS "Boost includes: " ${Boost_INCLUDE_DIRS} )
    message(STATUS "Boost libs: " ${Boost_LIBRARIES} )

#find_package (HDF5 COMPONENTS CXX)
#    message(STATUS "HDF5 libs: " ${HDF5_LIBRARIES} )

find_package (Eigen3 3.1)
    message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )



#find_package (PythonLibs REQUIRED)
#message("Include dirs of Python: " ${PYTHON_INCLUDE_DIRS} )
#message("Libs of Python: " ${PYTHON_LIBRARIES} )

include_directories(
    ${Boost_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_PATH}
    #${HDF5_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
    ./src
)

add_subdirectory(src)

option(Testing "Enable testing" ON)
if (Testing)
    add_subdirectory(tests)
    enable_testing()
endif (Testing)

# Build example
option(Examples "Build examples" OFF)
if (Examples)
    message(STATUS "Building examples")
    add_subdirectory(example)
endif (Examples)

# Install pkg-config file
configure_file(pomerol.pc.in pomerol.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pomerol.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)
