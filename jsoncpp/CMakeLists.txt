include_directories(./include)

set(json_sources src/json_reader.cpp src/json_writer.cpp src/json_value.cpp) 
add_library(pomerol_jsoncpp SHARED ${json_sources})

install(TARGETS pomerol_jsoncpp
    LIBRARY DESTINATION lib
    )
