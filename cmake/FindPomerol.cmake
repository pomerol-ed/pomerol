#  Try to find Pomerol. Once done this will define
#  Pomerol_FOUND - System has Pomerol
#  Pomerol_INCLUDE_DIRS - The Pomerol include directories
#  Pomerol_DEFINITIONS - Compiler switches required for using Pomerol

find_package(PkgConfig)
pkg_check_modules(PC_Pomerol QUIET pomerol)
set(Pomerol_DEFINITIONS ${PC_Pomerol_CFLAGS_OTHER})

#include
find_path(Pomerol_INCLUDE_DIR TwoParticleGF.h
          HINTS ${PC_Pomerol_INCLUDEDIR} ${PC_Pomerol_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} "${Pomerol_ROOT}/include/pomerol" 
         )
#lib
message (STATUS ${Pomerol_ROOT}/include/pomerol)
find_library(Pomerol_LIBRARY NAMES pomerol libpomerol
             HINTS ${PC_Pomerol_LIBDIR} ${PC_Pomerol_LIBRARY_DIRS} ${Pomerol_ROOT}/lib)

set(Pomerol_INCLUDE_DIRS ${Pomerol_INCLUDE_DIR} )
set(Pomerol_LIBRARIES ${Pomerol_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set Pomerol_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GFTools "No Pomerol found" Pomerol_LIBRARY Pomerol_INCLUDE_DIR)

mark_as_advanced(Pomerol_INCLUDE_DIR Pomerol_INCLUDE_DIRS Pomerol_LIBRARY Pomerol_LIBRARIES Pomerol_DEFINITIONS)
