
install( DIRECTORY .
         DESTINATION include
         FILES_MATCHING 
         PATTERN "*.h"
         PATTERN "*.hpp"
        )

configure_file(pomerol.h.in ${CMAKE_BINARY_DIR}/include/pomerol.h)
INSTALL(FILES ${CMAKE_BINARY_DIR}/include/pomerol.h DESTINATION include)


file(GLOB_RECURSE pomerol_headers RELATIVE ${CMAKE_CURRENT_LIST_DIR} "*.h" "*.hpp")
# generate the global include file
unset(PARSE_HEADERS)
foreach(include_file ${pomerol_headers})
    set(PARSE_HEADERS "${PARSE_HEADERS}\n#include <${include_file}>")
endforeach()

#message(STATUS $<JOIN:{PARSE_HEADERS}>)
#file(GENERATE OUTPUT "parse.h" CONTENT $<JOIN:${PARSE_HEADERS}>)

