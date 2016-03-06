#
# List of common macros for ${PROJECT_NAME} objects compilation
#

# Disable in-source builds
macro(no_source_builds)
    if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
        message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
    endif()
    # Print build directory
    message(STATUS "BUILD_DIR: ${CMAKE_BINARY_DIR}")
    # Print source directory
    message(STATUS "SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
    set(CMAKE_DISABLE_SOURCE_CHANGES ON)
    set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
endmacro(no_source_builds)

# Set global rpath
macro(fix_rpath)
    # RPATH fix
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
     set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
    else()
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif()

    #policy update CMP0042
    if(APPLE)
      set(CMAKE_MACOSX_RPATH ON)
    endif()
endmacro(fix_rpath)

# C++11
macro(set_cxx11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
        execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
        if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
            message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
        endif ()
    elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
        if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
            message( STATUS "Adding -stdlib=libc++ flag")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        endif()
    else ()
        message(FATAL_ERROR "Your C++ compiler does not support C++11.")
    endif ()
endmacro(set_cxx11)

# workaround definitions for different compilers
macro(compiler_workarounds)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
        set(BOLD_HYB_INTEL_BUGFIX TRUE)
        # disable "undefined template warning" frequently shown by intel
        add_definitions("-wd488")
    elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
        add_definitions("-Wno-sign-compare")
        add_definitions("-Wno-comment")
        set(BOLD_HYB_GNU_BUGFIX TRUE)
    endif()
endmacro(compiler_workarounds)

macro(set_linking_type)
option(BuildStatic "Build static libraries" OFF)
option(BuildShared "Build shared libraries" ON)
if (BuildStatic AND NOT BuildShared) 
    message(STATUS "Building static libraries")
    set(${PROJECT_CAPNAME}_BUILD_TYPE STATIC)
    set(BUILD_SHARED_LIBS OFF)
elseif(BuildShared AND NOT BuildStatic)
    message(STATUS "Building shared libraries")
    set(${PROJECT_CAPNAME}_BUILD_TYPE SHARED)
    set(BUILD_SHARED_LIBS ON)
else()
    message(FATAL_ERROR "Please choose EITHER BuildStatic OR BuildShared type of building libraries, NOT both")
endif()
endmacro(set_linking_type)

# 
# Dependencies
#

# ALPSCore - brings MPI, HDF5, boost 
macro(add_alpscore)
find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
    message(STATUS "ALPSCore includes:  ${ALPSCore_INCLUDES}")
    message(STATUS "ALPSCore libraries: ${ALPSCore_LIBRARIES}")
    target_include_directories(${PROJECT_NAME} PUBLIC ${ALPSCore_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${ALPSCore_LIBRARIES})
endmacro(add_alpscore)

# Eigen
macro(add_eigen3)
find_package (Eigen3 3.1 REQUIRED)
    message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )
    target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})
endmacro(add_eigen3)


# find and include gftools
macro(add_gftools)
find_package(GFTools QUIET)
    if (NOT GFTOOLS_FOUND)
        message(STATUS "Fetching gftools")
        include(ExtGFTools)
    endif()
    message(STATUS "GFTools includes : ${GFTOOLS_INCLUDE_DIR}") 
    include_directories(${GFTOOLS_INCLUDE_DIR})
endmacro(add_gftools)

# find, include and link to python
macro(add_python)
    find_package (PythonInterp)
    # tell cmake to find libs for found PythonInterp
    set(Python_ADDITIONAL_VERSIONS "${PYTHON_VERSION_STRING}")
    # There is a bug, involving cmake calling PYTHON framework on APPLE, which always happens to be the system python
    if (APPLE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_config_var(\"LIBDIR\")," OUTPUT_VARIABLE P_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_python_inc()," OUTPUT_VARIABLE P_INC OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(PYTHON_INCLUDE_DIR "${P_INC}")
        set(PYTHON_LIBRARY "${P_LIB}/libpython${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.dylib") # dirty hack
        unset(P_INC) 
        unset(P_LIB)
    endif(APPLE)

    find_package (PythonLibs REQUIRED)
    include_directories(${PYTHON_INCLUDE_DIRS})
    message(STATUS "Python includes: " ${PYTHON_INCLUDE_DIRS})
    message(STATUS "Python libraries: " ${PYTHON_LIBRARIES})
    #link_libraries(${PYTHON_LIBRARIES})
    set(OPENDF_HAS_PYTHON TRUE)
    list(APPEND boost_components python)
endmacro(add_python)

macro(add_fftw3)
find_package (fftw REQUIRED)
    message(STATUS "FFTW includes: " ${FFTW_INCLUDE_DIR} )
	target_include_directories(${PROJECT_NAME} PUBLIC ${FFTW_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${FFTW_LIBRARIES})
endmacro(add_fftw3)


# boost
macro(add_boost) # usage: add_boost(component1 component2...)
  find_package (Boost 1.54.0 COMPONENTS ${ARGV} REQUIRED)
  message(STATUS "Boost includes: ${Boost_INCLUDE_DIRS}" )
  message(STATUS "Boost libs: ${Boost_LIBRARIES}" )
  target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
endmacro(add_boost)

#
# Compilation/Installation shortcuts
#

# from alpscore
#macro(add_testing)
#  option(Testing "Enable testing" ON)
#  if (Testing)
#    enable_testing()
#    add_subdirectory(test)
#  endif (Testing)
#endmacro(add_testing)

# Compile tests in "test" subdir 
macro(add_testing)
option(Testing "Enable testing" ON)
if (Testing)
    include(EnableGtests) #defined in common/cmake
    enable_testing()
    add_subdirectory(test)
endif(Testing)
endmacro(add_testing)

# Build executables in "prog" subdirectory
macro(add_progs)
  foreach(src_ ${ARGV})
    set(prog_name "${src_}.${PROJECT_NAME}")
    add_executable(${prog_name} "prog/${src_}.cpp")
    target_link_libraries(${prog_name} ${PROJECT_NAME})
    set_target_properties(${prog_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
    install ( TARGETS ${prog_name} DESTINATION bin )
  endforeach()
endmacro(add_progs)

macro(add_mpi)
    find_package (MPI)
    message(STATUS "MPI includes: " ${MPI_CXX_INCLUDE_PATH} )
    message(STATUS "MPI CXX libs: " ${MPI_CXX_LIBRARIES} )
    message(STATUS "MPI definitions : " ${MPI_CXX_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
    add_definitions(${MPI_CXX_COMPILE_FLAGS}  ${MPI_C_COMPILE_FLAGS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES})
    target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_CXX_INCLUDE_PATH})
endmacro(add_mpi)

# Usage: add_this_package(srcs...)
# Compiles files in "src" into a library and installs it along with headers and cmake info
# The `srcs` are source file names in directory "src/"
# Defines ${PROJECT_NAME} target
# Exports alps::${PROJECT_NAME} target
function(add_this_package)
   # This is needed to compile tests:
   include_directories(
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_BINARY_DIR}/include
   )
  
  set(src_list_ "")
  foreach(src_ ${ARGV})
    list(APPEND src_list_ "src/${src_}.cpp")
  endforeach()
  add_library(${PROJECT_NAME} ${${PROJECT_CAPNAME}_BUILD_TYPE} ${src_list_})
  set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

  install(TARGETS ${PROJECT_NAME} 
          EXPORT ${PROJECT_NAME} 
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib
          INCLUDES DESTINATION include)
  install(EXPORT ${PROJECT_NAME} NAMESPACE ${PROJECT_NAME}:: DESTINATION share/${PROJECT_NAME})
  target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

  install(DIRECTORY include DESTINATION .
          FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx" PATTERN "*.h"
         )
endfunction(add_this_package)

# add interdependencies between ${PROJECT_NAME} modules
# Usage: add_opendf_module(pkgname1 pkgname2...)
# Sets variable ${PROJECT_NAME}_DEPENDS
macro(add_opendf_module)
    list(APPEND ${PROJECT_NAME}_DEPENDS ${ARGV})
    foreach(pkg_ ${ARGV})
        if (DEFINED BOLD_HYB_GLOBAL_BUILD)
            include_directories(${${pkg_}_INCLUDE_DIRS}) # this is needed to compile tests
            message(STATUS "${pkg_} includes: ${${pkg_}_INCLUDE_DIRS}" )
        else(DEFINED BOLD_HYB_GLOBAL_BUILD)
            string(REGEX REPLACE "^opendf_" "" pkgcomp_ ${pkg_})
            find_package(opendf QUIET COMPONENTS ${pkgcomp_})
            if (opendf_${pkgcomp_}_FOUND) 
              set(${pkg_}_LIBRARIES ${opendf_${pkgcomp_}_LIBRARIES})
            else()
              find_package(${pkg_} REQUIRED HINTS ${BOLD_HYB_ROOT})
            endif()
        endif (DEFINED BOLD_HYB_GLOBAL_BUILD)
        target_link_libraries(${PROJECT_NAME} PUBLIC ${${pkg_}_LIBRARIES})
        message(STATUS "${pkg_} libs: ${${pkg_}_LIBRARIES}")
    endforeach(pkg_)
endmacro(add_opendf_module)

# Generate documentation 
macro(gen_documentation)
  set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES} ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src" )
  option(Documentation "Build documentation" ON)
  if (Documentation)
    if (NOT BOLD_HYB_GLOBAL_BUILD)
        set(DOXYFILE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/include")
        set(DOXYFILE_IN "${PROJECT_SOURCE_DIR}/Doxyfile.in") 
        set(DOXYFILE_NAME ${PROJECT_NAME})
        set(DOXYFILE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc")
        include(UseDoxygen)
    else()
        set(DOXYFILE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/include")
        set(DOXYFILE_EXTRA_SOURCES "${DOXYFILE_EXTRA_SOURCES}")
        set(DOXYFILE_IN "${CMAKE_SOURCE_DIR}/Doxyfile.in") 
        set(DOXYFILE_NAME "${PROJECT_NAME}")
        set(DOXYFILE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc")
    endif()
  endif(Documentation)
endmacro(gen_documentation)

# Generate the automated configuration file with some of the compiler definitions
macro(gen_config_hpp)
  configure_file("${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/config.hpp.in" "${PROJECT_BINARY_DIR}/include/${PROJECT_NAME}/config.hpp")
  install(FILES "${PROJECT_BINARY_DIR}/include/${PROJECT_NAME}/config.hpp" DESTINATION include/${PROJECT_NAME}) 
endmacro(gen_config_hpp)

# Generation of pkg-config
macro(gen_pkg_config)
  # Generate pkg-config file
  configure_file("${PROJECT_SOURCE_DIR}/${PROJECT_NAME}.pc.in" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.pc")
  install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.pc" DESTINATION "lib/pkgconfig")
endmacro(gen_pkg_config)

# Function: generates package-specific CMake configs
# Arguments: [DEPENDS <list-of-dependencies>] [EXPORTS <list-of-exported-targets>]
# If no <list-of-dependencies> are present, the contents of ${PROJECT_NAME}_DEPENDS is used
# If no exported targets are present, alps::${PROJECT_NAME} is assumed.
function(gen_cfg_module)
    include(CMakeParseArguments) # arg parsing helper
    cmake_parse_arguments(gen_cfg_module "" "" "DEPENDS;EXPORTS" ${ARGV})
    if (gen_cfg_module_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Incorrect call of gen_cfg_module([DEPENDS ...] [EXPORTS ...]): ARGV=${ARGV}")
    endif()
    if (gen_cfg_module_DEPENDS)
        set(DEPENDS ${gen_cfg_module_DEPENDS})
    else()
        set(DEPENDS ${${PROJECT_NAME}_DEPENDS})
    endif()
    if (gen_cfg_module_EXPORTS)
        set(EXPORTS ${gen_cfg_module_EXPORTS})
    else()
        set(EXPORTS ${PROJECT_NAME}::${PROJECT_NAME})
    endif()
    configure_file("${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}ModuleConfig.cmake.in" 
                   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)
    configure_file("${PROJECT_SOURCE_DIR}/../cmake/${PROJECT_NAME}Config.cmake.in" 
                   "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)
    install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION "share/${PROJECT_NAME}/")
    install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION "share/${PROJECT_NAME}/")
    configure_file("${PROJECT_SOURCE_DIR}/../${PROJECT_NAME}.lmod.in" "${CMAKE_BINARY_DIR}/${PROJECT_NAME}.lmod")
endfunction()



