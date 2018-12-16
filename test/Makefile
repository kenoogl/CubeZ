# Makefile

target_arch=intel
PMdir=${HOME}/CubeZ/PMlib
PAPIdir= #/home/app/papi/5.5.1
fp_mode=float
#DEBUG_FLAG=-DDEBUG

ifeq ($(target_arch),intel)
  CXX = icpc
  #F90 = ifort
  #OPT = -xHOST -O3 -ipo -no-prec-div -fp-model fast=2 -qopenmp -qopt-report=3 -std=c++11
	OPT = -xHOST -O3 -ipo -no-prec-div -fp-model fast=2 -qopenmp -qopt-report=3 -D_NOAVX512
  CXXFLAGS = ${OPT} -I$(PMdir)/include -DDISABLE_MPI ${DEBUG_FLAG} #-I${PAPIdir}/include -DUSE_PAPI
  #F90FLAGS = ${OPT} -nofor-main -fpp -qopt-report=3
  #CXXFLAGS = -O3 -qopenmp -I$(PMdir)/PMlib/include -DDISABLE_MPI ${DEBUG_FLAG}
  #F90FLAGS = -O3 -nofor-main -qopenmp -fpp -qopt-report=3
  #LIBS = -lifport -lifcore -L$(PMdir)/lib -lPM
  LIBS = -L$(PMdir)/lib -lPM
  #LDFLAGS = $(PMdir)/lib/libpapi_ext.a ${PAPIdir}/lib/libpapi.a ${PAPIdir}/lib/libpfm.a
  LDFLAGS =
  ifeq ($(fp_mode),double)
    CXXFLAGS += -D_REAL_IS_DOUBLE_
    FFLAGS += -r8
  endif

### GNU
else ifeq ($(target_arch),gnu)
  CXX = g++
  F90 = gfortran
  CXXFLAGS = -O3 -fopenmp -I${PMdir}/PMlib/include -DDISABLE_MPI ${DEBUG_FLAG}
  F90FLAGS = -O3 -cpp -fopenmp -ffree-form
  LIBS = -L$(PMdir)/PMlib/lib -lPM
  LDFLAGS = -lgfortran
  ifeq ($(fp_mode),double)
    CXXFLAGS += -D_REAL_IS_DOUBLE_
    FFLAGS += -fdefault-real-8
  endif

### PGI
else ifeq ($(target_arch),pgi)
  CXX = pgc++
  F90 = pgfortran
  CXXFLAGS = -fastsse -mp -Mipa=fast,inline -Msafeptr --no_exceptions -Minfo=all -I${PMdir}/PMlib/include -DDISABLE_MPI ${DEBUG_FLAG}
  F90FLAGS = -fastsse -Minfo=all -Mipa=fast -Mbounds -Minfo=intensity -cpp -mp -Mfree
  LIBS = -L$(PMdir)/PMlib/lib -lPM
  LDFLAGS = -pgf90libs -fastsse -Mipa=fast,inline
  ifeq ($(fp_mode),double)
    CXXFLAGS += -D_REAL_IS_DOUBLE_
    FFLAGS += -r8
  endif

### FX
else ifeq ($(target_arch),fx)
  OPTFLAGS= -Kfast,ocl,preex,simd=2,array_private,parallel,openmp,uxsimd,optmsg=2 -V -Nsrc -x0 -Xg -Nfjcex -DDISABLE_MPI ${DEBUG_FLAG}
#-Kfast,parallel,ocl,preex,array_private,auto -Kopenmp,simd=2,uxsimd,optmsg=2 -V -Nsrc -Aw -Nobsfun -D_PM_WITHOUT_MPI_ -I${PMdir}/PMlib/include
#OPTFLAGS=-Kfast,parallel,ocl,preex,array_private,auto -Kopenmp,simd=2,uxsimd -Aw -Nobsfun -D__K_FAPP -D__K_FPCOLL
  CXX = FCCpx
  F90 = frtpx
  CXXFLAGS = $(OPTFLAGS)
  F90FLAGS = -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Cpp
  LIBS = -L$(PMdir)/PMlib/lib -lPM
  INCLUDES = -I${PMdir}/PMlib/include
  LDFLAGS = --linkfortran -Kopenmp
  ifeq ($(fp_mode),double)
    CXXFLAGS += -D_REAL_IS_DOUBLE_
    FFLAGS += -CcdRR8
  endif
endif

CSRC = test.cpp
FSRC =
SRCS = $(FSRC) $(CSRC)

.SUFFIXES: .o .cpp .f90

FOBJS = $(FSRC:.f90=.o)
COBJS = $(CSRC:.cpp=.o)
OBJS  = $(FOBJS) $(COBJS)

TARGET =  test

RM = rm

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(@) $(OBJS) $(LIBS) $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.f90.o:
	$(F90) $(F90FLAGS) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -Y *.cpp *.h
# DO NOT DELETE

test.o: test.cpp
