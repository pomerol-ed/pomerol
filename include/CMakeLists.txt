file(GLOB_RECURSE RELATIVE ${CMAKE_CURRENT_LIST_DIR} pomerol_headers "*.h")

install( FILES ${pomerol_headers}
         DESTINATION include)

#configure_file(pomerol.hpp.in ${CMAKE_BINARY_DIR}/include/pomerol.hpp)
