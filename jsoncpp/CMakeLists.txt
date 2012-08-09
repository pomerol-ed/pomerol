include_directories(./include)

set(json_sources src/json_reader.cpp src/json_writer.cpp src/json_value.cpp) 
add_library(jsoncpp SHARED ${json_sources})
