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


add_executable("test1" "test1.cpp")
target_link_libraries("test1" CZ FCORE)

add_executable("test2" "test2.cpp")
target_link_libraries("test2" CZ FCORE)

add_executable("test3" "test3.cpp")
target_link_libraries("test3" CZ FCORE)

###
