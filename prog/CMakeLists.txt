include_directories (. ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    pomerol
)

add_executable(pomerolDiag DevelopersExample.cpp)
install ( TARGETS pomerolDiag DESTINATION bin )
