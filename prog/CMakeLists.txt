include_directories (. ../jsoncpp/include ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    jsoncpp
    pomerol
)

add_executable(pomerolDiag DevelopersExample.cpp)
install ( TARGETS pomerolDiag DESTINATION bin )
