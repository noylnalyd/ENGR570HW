cmake_minimum_required(VERSION 3.10)

# set the project name
project(HW2 Fortran)

# add the executable
SET(CMAKE_FC_COMPILER /use/bin/gfortran)
add_executable(HW2 HW2.f90 mm_io.f90)