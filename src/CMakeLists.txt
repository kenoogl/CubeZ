###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################

add_subdirectory(cz_cpp)
add_subdirectory(cz_f90)

include_directories(
       ${PROJECT_BINARY_DIR}/src/cz_cpp  # czVersion.h
       ${PROJECT_SOURCE_DIR}/src/cz_cpp
       ${PROJECT_SOURCE_DIR}/src/cz_f90
)

link_directories(
       ${PROJECT_BINARY_DIR}/src/cz_cpp
       ${PROJECT_BINARY_DIR}/src/cz_f90
)


IF(with_PM)
  include_directories(${PM_INC})
  link_directories(${PM_LIB})
endif()


IF(with_PAPI)
  include_directories(${PAPI_INC})
  link_directories(${PAPI_LIB})
endif()

if(with_CBR)
  include_directories(${CBR_INC})
  link_directories(${CBR_LIB})
endif()


set(cz_src main.cpp)

if(with_MPI)
  set(cz_target "cz-mpi")
else()
  set(cz_target "cz")
endif()



add_executable(${cz_target} ${cz_src})

if(with_ACC)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set_target_properties(${cz_target} PROPERTIES LINKER_LANGUAGE Fortran)
  endif()
else()
  set_target_properties(${cz_target} PROPERTIES LINKER_LANGUAGE CXX)
endif()

target_link_libraries(${cz_target} CZ FCORE)


install(TARGETS ${cz_target} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


###
