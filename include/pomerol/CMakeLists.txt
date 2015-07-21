file(GLOB pomerol_headers "*.h")

install( FILES ${pomerol_headers}
         DESTINATION include/pomerol )

#configure_file(pomerol.hpp.in ${CMAKE_BINARY_DIR}/include/pomerol.hpp)
