set(default_progs
    anderson
    hubbard2d
    quantum_model)

# Find or fetch gftools
add_gftools()

message(STATUS "Building executables")
if(NOT progs_list)
    set(progs_list ${default_progs})
endif()

# TODO
#foreach(model ${progs_list})
#  list (FIND default_progs ${model} _index)
#  if (${_index} GREATER -1)
#    string(TOUPPER ${model} def)
#    set(def "POMEROL_${def}")
#    add_prog(${model} ${def})
#  else()
#    message(WARNING "Unknown program \"${model}\" has been skipped.")
#  endif()
#endforeach()
