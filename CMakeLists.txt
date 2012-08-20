cmake_minimum_required (VERSION 2.6)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project (pomerol)

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

find_package (Boost COMPONENTS python)
find_package (HDF5 COMPONENTS CXX)
find_package (Eigen3 3.1)

message(STATUS "Boost includes: " ${Boost_INCLUDE_DIRS} )
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )
message(STATUS "Boost libs: " ${Boost_LIBRARIES} )
message(STATUS "HDF5 libs: " ${HDF5_LIBRARIES} )

#find_package (PythonLibs REQUIRED)
#message("Include dirs of Python: " ${PYTHON_INCLUDE_DIRS} )
#message("Libs of Python: " ${PYTHON_LIBRARIES} )


include_directories(
    ${Boost_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
    ./jsoncpp/include
    ./src
)

add_subdirectory(jsoncpp)
add_subdirectory(src)
add_subdirectory(tests)
enable_testing()


