# CubeZ

## OUTLINE

CubeZ is a platform for testing iterative solvers.


## Copyright
- Copyright (c) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University. All rights reserved.



## Prerequisite

- Cmake
- MPI library (if parallel)
- PMlib
- PAPI (Optional)
- CBrick (In case of parallel)


## How to build

### Build

~~~
$ export HOME=/hogehoge
$ mkdir build
$ cd build
$ cmake [options] ..
$ make
$ sudo make install
~~~


### Options

`-D INSTALL_DIR=` *Install_directory*
>  Specify the directory that this library will be installed. Built library is
   installed at `install_directory/lib` and the header files are placed at
   `install_directory/include`.

`-D enable_OPENMP=` {yes | no}
>  This option makes OpenMP directives effect. Default is yes.

`-D with_MPI=` {yes | no}
>  If you use an MPI library, specify `with_MPI=yes` (default).

`-D real_type=` {float | double}
>  Specify the type of floating point. If this option is omitted, the default is float.

`-D with_PM=` {*Installed_Directory* | OFF}
> Specify the directory path that PMlib is installed, or OFF.

`-D with_FAPP=` {OFF | ON}
> Specify using Fujitsu profiler.

`-D with_CBR=` *Installed_Directory*
> Specify the directory path that CBrick is installed.

`-D with_PAPI=` *Installed_Directory*
> Specify the directory path that PAPI is installed.

`-D with_SIMD=` {OFF |256|512}
> Specify SIMD length. The default is OFF. If you want to use AVX512 specify 512.

`-D with_ACC=`{off|Pascal|Volta}
> Specify using open acc directives.

`-D with_Ftrace=` (off | on)
> In the case of Aurora, if you want to use Ftrace option, specify turn on this option.

`-D enable_VectorReduction=` (on | off)
> Vectorizable reduction expression. Default is yes. This is effective for Aurora.


### Default settng
~~~
with_MPI = ON
enable_OPENMP = ON
real_type = float
with_PAPI = OFF
with_FAPP = OFF
with_SIMD = OFF
with_Ftrace = OFF
with_ACC = OFF
enable_VectorReduction = ON
~~~


## Configure Examples

`$ export CZ_HOME=hogehoge`

In following examples, assuming that TextParser, PMlib, and CBrick are installed under the CZ_HOME directory. If not, please specify applicable directory paths.


### Serial build

When even compiling without MPI, specify `-D with_CBR=OFF` to suppress linking to CBrick library.
In the case of some Intel compiler environment, please specify environment variables before compilation.
`export CC=icc CXX=icpc F90=ifort FC=ifort`


### Mac

#### gnu serial /w PMlib /wo PAPI

~~~
$ module load gcc
$ echo $CC $CXX $FC $F90
gcc-9 g++-9 gfortran-9 gfortran-9

$ cmake -DINSTALL_DIR=${HOME}/CubeZ/CZ \
-Denable_OPENMP=ON \
-Dwith_PM=${HOME}/CubeZ/PMlib \
-Dwith_SIMD=256 ..
~~~

#### PGI serial /w PMlib /wo PAPI

~~~
$ module load pgi/19.4
$ export CC=pgcc CXX=pgc++ F90=pgf90 FC=pgf90

$ cmake -DINSTALL_DIR=${HOME}/CubeZ/CZ \
-Denable_OPENMP=ON \
-Dwith_PM=${HOME}/CubeZ/PMlib_PGI \
-Dwith_SIMD=256 ..
~~~


#### INTEL serial /w PMlib /wo PAPI

~~~
$ module load intel
$ export CC=icc CXX=icpc F90=ifort FC=ifort

$ cmake -DINSTALL_DIR=${HOME}/CubeZ/CZ \
-Denable_OPENMP=ON \
-Dwith_PM=${HOME}/CubeZ/PMlib \
-Dwith_SIMD=256 ..
~~~



### ITO subsystem A/B

####  Intel serial /w PMlib /w PAPI

~~~
$ module load intel/2018
$ export CC=icc CXX=icpc F90=ifort FC=ifort

$ cmake -DINSTALL_DIR=${HOME}/CZ \
-Denable_OPENMP=ON \
-Dwith_PM=${HOME}/opt/PMlib/intel-2018_papi-gcc-4.8.5 \
-Dwith_PAPI=${HOME}/opt/PAPI/gcc-4.8.5 \
-Dwith_SIMD=256 ..
~~~

#### Intel MPI /w PMlib /w PAPI

- 2019-12-25 性能はとれるが、計算結果が不正で反復回数が増加する

~~~
$ module load intel/2018.3 openmpi/3.1.3-nocuda-intel18.3

>> PAPI-5.7.0
$ export CC=mpiicc CXX=mpiicpc F90=mpiifort FC=mpiifort
$ ./configure --prefix=${HOME}/CZ_INTL/PAPI
$ make
$ make install

>> PMlib-6.4.5
$ cmake -DINSTALL_DIR=${HOME}/CZ_INTL/PMlib \
-Denable_OPENMP=ON \
-Dwith_MPI=ON \
-Dwith_PAPI=${HOME}/CZ_INTL/PAPI \
-Denable_Fortran=no \
-Dwith_example=no \
-Dwith_OTF=no \
-Denable_PreciseTimer=yes ..

>> CBrick-1.4.3
$ cmake -DINSTALL_DIR=${HOME}/CZ_INTL/CBrick \
      -Denable_OPENMP=ON ..

>> CubeZ-1.2.9
$ cmake -DINSTALL_DIR=${HOME}/CZ_INTL/CZ \
-Denable_OPENMP=ON \
-Dwith_MPI=ON \
-Dwith_PM=${HOME}/CZ_INTL/PMlib \
-Dwith_PAPI=${HOME}/CZ_INTL/PAPI \
-Dwith_SIMD=256 \
-Dwith_CBR=${HOME}/CZ_INTL/CBrick ..
~~~

#### Fujitsu TCS serial /w PMlib /w PAPI /w FAPP

- 2020-01-07

~~~
>> PAPI-5.7.0
$ module load gcc
$ export CC=gcc F77=gfortran
$ ./configure --prefix=${HOME}/CZ_TCS/PAPI
$ make
$ make install

$ export CC=fcc CXX=FCC F90=frt FC=frt

>> PMlib-6.4.5
$ cmake -DINSTALL_DIR=${HOME}/CZ_TCS/PMlib \
-DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_ITO_TCS.cmake \
-Denable_OPENMP=yes \
-Dwith_MPI=no \
-Denable_Fortran=no \
-Dwith_example=no \
-Dwith_PAPI=${HOME}/CZ_TCS/PAPI \
-Dwith_OTF=no \
-Denable_PreciseTimer=OFF ..

>> CubeZ-1.2.9
$ cmake -DINSTALL_DIR=${HOME}/CZ_TCS/CZ \
-DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_ITO_TCS.cmake \
-Denable_OPENMP=ON \
-Dwith_PAPI=${HOME}/CZ_TCS/PAPI \
-Dwith_PM=${HOME}/CZ_TCS/PMlib \
-Dwith_FAPP=ON \
-Dwith_SIMD=256 ..
~~~

#### Fujitsu TCS parallel /wo PMlib /w FAPP

- 2020-01-07

~~~
>> PAPI-5.7.0
$ module load gcc
$ export CC=gcc F77=gfortran
$ ./configure --prefix=${HOME}/CZ_TCS/PAPI
$ make
$ make install

$ export CC=mpifcc CXX=mpiFCC F90=mpifrt FC=mpifrt

>> PMlib-6.4.5
$ cmake -DINSTALL_DIR=${HOME}/CZ_TCS/PMlib \
-DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_ITO_TCS.cmake \
-Denable_OPENMP=ON \
-Dwith_MPI=ON \
-Denable_Fortran=OFF \
-Dwith_example=OFF \
-Dwith_PAPI=${HOME}/CZ_TCS/PAPI \
-Dwith_OTF=OFF \
-Denable_PreciseTimer=OFF ..

>> CBrick-1.4.3
$ cmake -DINSTALL_DIR=${HOME}/CZ_TCS/CBrick \
-DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_ITO_TCS.cmake \
-Denable_OPENMP=ON ..

>> CubeZ-1.2.9
$ cmake -DINSTALL_DIR=${HOME}/CZ_TCS/CZ \
-DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_ITO_TCS.cmake \
-Denable_OPENMP=ON \
-Dwith_MPI=ON \
-Dwith_PAPI=${HOME}/CZ_TCS/PAPI \
-Dwith_PM=${HOME}/CZ_TCS/PMlib \
-Dwith_FAPP=ON \
-Dwith_SIMD=256 \
-Dwith_CBR=${HOME}/CZ_TCS/CBrick ..
~~~




#### PGI OpenACC with PAPI

~~~
$ module load pgi/19.4
$ export CC=pgcc CXX=pgc++ F90=pgf90 FC=pgf90
~~~

~~~
$ cmake -DINSTALL_DIR=${HOME}/CZ \
-Dwith_MPI=no \
-Dreal_type=float \
-Denable_OPENMP=no \
-Dwith_PM=${HOME}/opt/PMlib/pgi-19.4_papi-gcc-4.8.5 \
-Dwith_SIMD=256 \
-Dwith_PAPI=${HOME}/opt/PAPI/gcc-4.8.5 \
-Dwith_ACC=Pascal \
-Dwith_CBR=OFF ..
~~~

#### PGI CPU OpenMP with PAPI

~~~
$ cmake -DINSTALL_DIR=${HOME}/CZ \
-Dwith_MPI=no \
-Dreal_type=float \
-Denable_OPENMP=yes \
-Dwith_PM=${HOME}/opt/PMlib/pgi-19.4_papi-gcc-4.8.5 \
-Dwith_SIMD=256 \
-Dwith_PAPI=${HOME}/opt/PAPI/gcc-4.8.5 \
-Dwith_ACC=OFF \
-Dwith_CBR=OFF ..
~~~




### Aurora without PAPI

~~~
$ export CC=ncc CXX=nc++ F90=nfort FC=nfort
~~~

~~~
$ cmake -DINSTALL_DIR=${HOME}/CZ \
-DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_NEC_Aurora.cmake \
-Dwith_MPI=OFF \
-Dwith_PM=OFF \
-Dwith_CBR=OFF \
-Dwith_Ftrace=ON -Denable_VectorReduction=OFF ..
~~~



##### Note
- On Fujitsu machines(K, fx100), confirm appropriate directory path for compiler environment.
- Before building, execute following command to clean for sure. `$ make distclean`


## Building PAPI by gcc
  
  ~~~
  $ gcc --version
  gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-11)
  
  $ export CC=gcc F77=gfortran
  $ cd papi-5.7.0/src
  $ ./configure --prefix=${HOME}/opt/PAPI/gcc-4.8.5
  $ make
  $ make install
  ~~~


## Building PMlib by PGI
  gccでビルドしたPAPIオブジェクトをリンクする
  
  ~~~
  $ cd PMlib-6.4.2
  $ mkdir build
  $ cd build
  $ module load pgi/19.4
  $ export CC=pgcc CXX=pgc++ F90=pgf90 FC=pgf90
  
  $ cmake -DINSTALL_DIR=${HOME}/opt/PMlib/OMP-pgi-19.4_papi-gcc-4.8.5 \
	-Denable_OPENMP=yes \
	-Dwith_MPI=no \
	-Denable_Fortran=no \
	-Dwith_example=no \
	-Dwith_PAPI=${HOME}/opt/PAPI/gcc-4.8.5 \
	-Dwith_OTF=no \
    -Denable_PreciseTimer=yes ..
  $ make
  $ make install
  ~~~


## Contributors

- Kenji Ono


### Comment of precision of Fortran
For example, the `-CcdRR8` option for fortran preprocessor convert variables, functions, and constants higher precision version in source code. Thus, functions in source code is described by floating version.

## 使い方

~~~
$ ./cz gsz_x, gsz_y, gsz_z, linear_solver, IterationMax coef [gdv_x, gdv_y, gdv_z]
$ ./cz 124 124 124 sor2sma 10000 1.5
$ ./cz 124 124 124 sor2sma_maf 10000 1.5
$ ./cz 124 124 124 pbicgstab 10000 1.5 {jacobi, psor, sor2sma, pcr}
$ ./cz 124 124 124 pbicgstab_maf 10000 1.5 {jacobi, psor, sor2sma, pcr}
$ ./cz 124 124 124 pcr 10000 1.5
$ ./cz 124 124 124 pcr_maf 10000 1.5
~~~
 - gsz_x, gsz_y, gsz_z  全計算領域の要素数
 - linear_solver        線形ソルバの指定
   - jacobi
   - sor
   - sor2sma
   - pbicgstab
   - pcr
 - IterationMax         最大反復回数
 - coef  緩和/加速係数
 - gdv_x, gdv_y, gdv_z  領域分割数の指定、指定しない場合には自動分割
