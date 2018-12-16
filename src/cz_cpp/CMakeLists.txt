###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################


include_directories(
       ${PROJECT_BINARY_DIR}/src/cz_cpp  # czVersion.h
       ${PROJECT_SOURCE_DIR}/src/cz_cpp
       ${PL_INC}
)

link_directories(
      ${PL_LIB}
)

IF(with_PM)
  include_directories(${PM_INC})
  link_directories(${PM_LIB})
endif()

IF(with_PAPI)
  include_directories(${PAPI_INC})
  link_directories(${PAPI_LIB})
  # link static library
  set(PAPI_STATIC_LIB ${PAPI_LIB}/libpapi.a
                      ${PAPI_LIB}/libpfm.a
                      )
endif()

if(with_CBR)
  include_directories(${CBR_INC})
  link_directories(${CBR_LIB})
endif()



set(cz_files
       cz_Evaluate.cpp
       cz_miscel.cpp
       cz_Poisson.cpp
       cz_comm.cpp
       cz_dsolver.cpp
       cz_simd.cpp
)

add_library(CZ STATIC ${cz_files})

if(with_MPI)
  if(with_CBR)
    if(real_type STREQUAL "double")
      set(cz_libs "-lCBrick_d")
    else()
      set(cz_libs "-lCBrick_f")
    endif()
  endif()

  if(with_PM)
    set(cz_libs "${cz_libs} -lPMmpi")
  endif()

else()
  if(with_PM)
    set(cz_libs "-lPM")
  endif()
endif()

if(with_PAPI)
  if(TARGET_ARCH STREQUAL "FX100")
  # on AICS FX100 static PAPI libs does not match the other system libs.
  # so the following dynamic libs must be used. bonehead Fujitsu...
    set(cz_libs "${cz_libs} -lpapi_ext -lpapi -lpfm")
    ## target_link_libraries(example1 -lpapi_ext -Wl,'-lpapi,-lpfm')
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  ## GNU linker ld options are better passed using -Wl, rather than -B compiler phase option
    #set(cz_libs "${cz_libs} -lpapi_ext -lpapi -lpfm")
    set(cz_libs "${cz_libs} -lpapi_ext ${PAPI_STATIC_LIB}")
    ## target_link_libraries(${cz_target} -lpapi_ext -Wl,'-Bstatic,-lpapi,-lpfm,-Bdynamic')
  else()
    ## Intel compiler & others
    # set(cz_libs "${cz_libs} -lpapi_ext -lpapi -lpfm")
    set(cz_libs "${cz_libs} -lpapi_ext ${PAPI_STATIC_LIB}")
  endif()
endif()

target_link_libraries(CZ "${cz_libs}")

###
