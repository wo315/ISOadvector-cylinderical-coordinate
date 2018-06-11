# CMake project file for NORMA

cmake_minimum_required (VERSION 2.6)
project (NORMA)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran*")
#  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fcray-pointer -fno-f2c -O3 -mcmodel=medium")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -fcray-pointer -mcmodel=medium")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-mcmodel=medium -f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran*")


# build executables
# set (NMPROGRAMS "diagstd" "diagrtb" "proj_modes_bin" "pdbmat")
# set (EXECUTABLES "NORMA.exe" ${NMPROGRAMS})
# set (SCRIPTS "gen_pert.sh" "pert_multi_mode.sh")

# add_executable ("NORMA.exe" "NORMA.f")
# foreach (p ${NMPROGRAMS})
#   add_executable (${p} "${p}.f")
# endforeach (p)

# install executables and scripts
# install (TARGETS ${EXECUTABLES}
#          RUNTIME DESTINATION "bin")
# install (PROGRAMS ${SCRIPTS}
#          DESTINATION "bin")

set(SRCS
  main.f90 
  par.f90
  set.f90
  initial.f90 
  ../dagmg.f90 
  ../acceleration.f90 
  ../mumps.f90 
  ../ADS.f90 
  ../next.f90 
  ../NS_equation.f90 
  ../blas.f par.f90 
  ../boundary.f90 
  ../property.f90 
  ../solver.f90 
  ../lapack.f 
  ../isoadvector.f90
  ../isolib.f90
  ../mesh.f90 
  ../HRscheme.f90 
  ../check.f90 
  ../CSF.f90
  ../output.f90
  )
add_executable(main ${SRCS})

target_link_libraries(main stdc++)
