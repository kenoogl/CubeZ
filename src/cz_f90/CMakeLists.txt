###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################

include_directories(
       ${CBR_INC}
)

set(cz_files
       cz_solver.f90
       cz_blas.f90
       cz_utility.f90
       cz_lsor.f90
       obsolete.f90
       cz_maf.f90
)

add_library(FCORE STATIC ${cz_files})

if(with_ACC)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    SET_TARGET_PROPERTIES(FCORE PROPERTIES LINKER_LANGUAGE Fortran)
  endif()
else()
  SET_TARGET_PROPERTIES(FCORE PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(FCORE ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
endif()

###
