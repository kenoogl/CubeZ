# CubeZ

## OUTLINE

CubeZ is a platform for testing iterative solvers.


## Copyright
- Copyright (c) 2018 Research Institute for Information Technology(RIIT), Kyushu University. All rights reserved.



## Prerequisite

- Cmake
- MPI library (if parallel)
- PMlib
- CBrick (Even compiling without MPI, a header file CB_Define.h is essential. See serial build.)


## How to build

### Build

~~~
$ export CZ_HOME=/hogehoge
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

`-D with_CBR=` *Installed_Directory*

> Specify the directory path that CBrick is installed.

`-D with_PAPI=` *Installed_Directory*

> Specify the directory path that PAPI is installed.


## Configure Examples

`$ export CZ_HOME=hogehoge`

In following examples, assuming that TextParser, PMlib, and CBrick are installed under the CZ_HOME directory. If not, please specify applicable directory paths.

### INTEL/GNU compiler

~~~
$ cmake -DINSTALL_DIR=${CZ_HOME}/CubeZ \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_PM=${CZ_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_CBR=${CZ_HOME}/CBrick ..
~~~

#### Note
In case of some Intel compiler environment, please specify environment variables before compilation.
`export CC=icc CXX=icpc F90=ifort FC=ifort`


### FUJITSU compiler / FX100, K computer on login nodes (Cross compilation) and Fujitsu TCS environment for intel PC

~~~
$ cmake -DINSTALL_DIR=${CZ_HOME}/RAinWATER \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx100.cmake \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_PM=${CZ_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_CBR=${CZ_HOME}/CBrick ..

$ cmake -DINSTALL_DIR=${CZ_HOME}/RAinWATER \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_PM=${CZ_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_CBR=${CZ_HOME}/CBrick ..

$ cmake -DINSTALL_DIR=${CZ_HOME}/RAinWATER \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_intel_F_TCS.cmake \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_PM=${CZ_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_CBR=${CZ_HOME}/CBrick ..
~~~

##### Note
- On Fujitsu machines(K, fx100), confirm appropriate directory path for compiler environment.
- Before building, execute following command to clean for sure. `$ make distclean`


### Serial build
When even compiling without MPI, the hearer file CB_Define.h in CBrick library is required.
Copy CB_Define.h form CBrick directory into cz_cpp, then make. In this case, specify `-D with_CBR=OFF`, then CBrick library is not linked.


## Contributors

- Kenji Ono


### Comment of precision of Fortran
For example, the `-CcdRR8` option for fortran preprocessor convert variables, functions, and constants higher precision version in source code. Thus, functions in source code is described by floating version.

## 使い方

~~~
$ ./cz gsz_x, gsz_y, gsz_z, linear_solver, IterationMax [gdv_x, gdv_y, gdv_z]
~~~
 - gsz_x, gsz_y, gsz_z  全計算領域の要素数
 - linear_solver        線形ソルバの指定
   - jacobi
   - sor
   - sor2sma
   - pbicgstab
   - lsor
 - IterationMax         最大反復回数
 - gdv_x, gdv_y, gdv_z  領域分割数の指定、指定しない場合には自動分割
