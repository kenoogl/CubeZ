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

set(tgt1 "test1")
add_executable(${tgt1} "${tgt1}.cpp")
target_link_libraries(${tgt1} CZ FCORE)

set(tgt2 "test2")
add_executable(${tgt2} "${tgt2}.cpp")
target_link_libraries(${tgt2} CZ FCORE)

set(tgt3 "test3")
add_executable(${tgt3} "${tgt3}.cpp")
target_link_libraries(${tgt3} CZ FCORE)

set(tgt4 "test4")
add_executable(${tgt4} "${tgt4}.cpp")
target_link_libraries(${tgt4} CZ FCORE)

set(tgt5 "test5")
add_executable(${tgt5} "${tgt5}.cpp")
target_link_libraries(${tgt5} CZ FCORE)

###
