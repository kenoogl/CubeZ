###################################################################################
#
# CubeZ
#
# Copyright (C) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################

##
## Compile option selector
##


macro (AddOptimizeOption)
  # from https://github.com/SX-Aurora/CMake-toolchain-file
  if(TARGET_ARCH STREQUAL "NEC_Aurora_VE")
    set(CMAKE_Fortran_COMPILER /opt/nec/ve/bin/nfort CACHE FILEPATH "Aurora Fortran compiler")
    set(CMAKE_CXX_COMPILER /opt/nec/ve/bin/nc++ CACHE FILEPATH "Aurora C++ compiler")
    set(CMAKE_C_COMPILER /opt/nec/ve/bin/ncc CACHE FILEPATH "Aurora C compiler")
    set(CMAKE_LINKER /opt/nec/ve/bin/nld CACHE FILEPATH "Aurora linker")
    set(CMAKE_AR /opt/nec/ve/bin/nar CACHE FILEPATH "Aurora archiver")
    set(CMAKE_RANLIB /opt/nec/ve/bin/nranlib CACHE FILEPATH "Aurora ranlib")
    set(CMAKE_CXX_FLAGS "-O3 -proginf")
    set(CMAKE_Fortran_FLAGS "-fpp -Wall -O3 -proginf -report-all -fdiag-parallel=2 -fdiag-vector=2 -std=f95 -cxxlib -static-nec -mretain-none")
    # In order to link objects by Fortran driver, add '-cxxlib' to CMAKE_Fortran_FLAGS

  elseif (TARGET_ARCH STREQUAL "ITO_TCS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,parallel,optmsg=2 -V -Xg")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,parallel,optmsg=2 -Xg")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Cpp -Kfast,parallel,optmsg=2 -V")

  elseif (USE_F_TCS STREQUAL "YES")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg -Nfjcex")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Cpp -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Qt")
    # -Xg   : gcc compatible flag
    # -fPIC : PIC flag

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")
    set(CMAKE_Fortran_FLAGS "-O3 -Wall -cpp")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopt-report=5 -std=c++11 -restrict  -DMPICH_IGNORE_CXX_SEEK")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -qopt-report=5")
    set(CMAKE_Fortran_FLAGS "-O3 -qopt-report=5 -inline-forceinline -fpp")

    # optimization for Intel AVX*
    if (with_SIMD STREQUAL "256")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX2 -D_SIMD_256")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xCORE-AVX2")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xCORE-AVX2 -align array32byte")
    elseif (with_SIMD STREQUAL "512")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512 -D_SIMD_512 -qopt-zmm-usage=high")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xCORE-AVX512 -D_SIMD_512 -qopt-zmm-usage=high")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xCORE-AVX512 -align array64byte -qopt-zmm-usage=high")
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -Mipa=fast,inline -O3 -g -Minfo=intensity,vect,mp,acc -cpp=mm -tp=skylake")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fast -Mipa=fast,inline -O3 -g -Minfo=intensity,mp,vect,acc -cpp=mm -tp=skylake")
    set(CMAKE_Fortran_FLAGS "-fast -Mipa=fast,inline -O3 -g -Minfo=intensity,mp,vect,acc -cpp=mm -tp=skylake")

    # optimization for Intel AVX*
    if (with_SIMD STREQUAL "256")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mvect=simd:256 -D_SIMD_256")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Mvect=simd:256")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mvect=simd:256")
    elseif (with_SIMD STREQUAL "512")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mvect=simd:512 -D_SIMD_512")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Mvect=simd:512")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mvect=simd:512")
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -cpp -Mfree")

  else()
    message("using default option")
  endif()
endmacro()


macro (FreeForm)
  if(TARGET_ARCH STREQUAL "NEC_Aurora_VE")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form")

  elseif(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  elseif (USE_F_TCS STREQUAL "YES")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Free")

  elseif(TARGET_ARCH STREQUAL "ITO_TCS")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Free")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfree")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfree")

  endif()
endmacro()


macro(C99)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endmacro()


macro(CPP11)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif()
endmacro()


macro(checkOpenMP)
  if(enable_OPENMP)
    if(TARGET_ARCH STREQUAL "NEC_Aurora_VE")
     set(OpenMP_Fortran_FLAGS "-fopenmp" CACHE STRING "Flag to enable OpenMP")
     set(OpenMP_CXX_FLAGS "-fopenmp" CACHE STRING "Flag to enable OpenMP")
     set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "Flag to enable OpenMP")

    elseif(USE_F_TCS STREQUAL "YES")
      set(OpenMP_C_FLAGS "-Kopenmp")
      set(OpenMP_CXX_FLAGS "-Kopenmp")
      set(OpenMP_Fortran_FLAGS "-Kopenmp")

    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(OpenMP_C_FLAGS "-fopenmp")
      set(OpenMP_CXX_FLAGS "-fopenmp")
      set(OpenMP_Fortran_FLAGS "-fopenmp")

    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      set(OpenMP_C_FLAGS "-qopenmp")
      set(OpenMP_CXX_FLAGS "-qopenmp")
      set(OpenMP_Fortran_FLAGS "-qopenmp")

    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
      set(OpenMP_C_FLAGS "-mp")
      set(OpenMP_CXX_FLAGS "-mp")
      set(OpenMP_Fortran_FLAGS "-mp")

    else()
      find_package(OpenMP REQUIRED)
    endif()

    # OpenMP_*_FLAGSにはfind_package(OpenMP REQUIRED)でオプションフラグが設定される
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endmacro()


macro(precision)
  if(real_type STREQUAL "OFF")
  # nothing, default is float
  set(real_type "float")

  elseif(real_type STREQUAL "float")
  # nothing

  elseif(real_type STREQUAL "double")
    ADD_DEFINITIONS(-D_REAL_IS_DOUBLE_)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_REAL_IS_DOUBLE_")

    if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -CcdRR8")

    elseif(CMAKE_Fortran_COMPILER MATCHES ".*frt$")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -CcdRR8")

    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")

    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")

    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
    endif()

  else() # neither 'float' nor 'double'
    message("@@@@@@@@@@@")
    message("FATAL ERROR : Invalid floating type : ${real_type}")
    message("@@@@@@@@@@@")
  ENDIF()
endmacro()
