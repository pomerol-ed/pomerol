# Install header files
install(DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.hpp")

# Prepare and install pomerol/Version.hpp
if(USE_OPENMP AND OPENMP_FOUND)
    set(POMEROL_USE_OPENMP ON)
endif()
configure_file("pomerol/Version.hpp.in" "pomerol/Version.hpp")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pomerol/Version.hpp"
        DESTINATION include/pomerol)
