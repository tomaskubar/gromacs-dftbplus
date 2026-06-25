#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# \author Victor Holanda <victor.holanda@cscs.ch>
# \author Joe Jordan <ejjordan@kth.se>
# \author Prashanth Kanduri <kanduri@cscs.ch>
# \author Sebastian Keller <keller@cscs.ch>

set(argon "argon-forces-integration")
add_executable(${argon} argon-forces-integration.cpp)

set(methane "methane-water-integration")
add_executable(${methane} methane-water-integration.cpp)

if (GMX_QMMM_DFTBPLUS)
    LINK_DIRECTORIES(GMX_QMMM_DFTBPLUS_LIB)
    add_library(libdftbplus STATIC IMPORTED)
    set_target_properties(libdftbplus PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libdftbplus.a)
#   add_library(libdftd3 STATIC IMPORTED)
#   set_target_properties(libdftd3 PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libdftd3.a)
#   add_library(libmudpack STATIC IMPORTED)
#   set_target_properties(libmudpack PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libmudpack.a)

    # search for the gfortran and gomp libraries in the system
    # first look for possible library locations
    unset(_lib_env)
    string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
    list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    list(REMOVE_DUPLICATES _lib_env)
    set(PATH_TO_LOOK_FOR "${_lib_env}")
    # try to find the gfortran lib in the given paths and add it
    find_library(GFORTRAN
        NAMES libgfortran.so.5 libgfortran.so.4
        HINTS ${PATH_TO_LOOK_FOR})
    add_library(gfortran SHARED IMPORTED)
    set_target_properties(gfortran PROPERTIES IMPORTED_LOCATION ${GFORTRAN})
    # try to find the gomp lib in the given paths and add it
    find_library(GOMP
        NAMES libgomp.so.1 libgomp.so
        HINTS ${PATH_TO_LOOK_FOR})
    add_library(gomp SHARED IMPORTED)
    set_target_properties(gomp PROPERTIES IMPORTED_LOCATION ${GOMP})

    set_target_properties(${argon} PROPERTIES LINK_FLAGS
      "-Wl,--whole-archive ${GMX_QMMM_DFTBPLUS_LIB}/libdftbplus.a -Wl,--no-whole-archive")
    target_link_libraries(${argon} PRIVATE nblib libdftbplus ${LINEAR_ALGEBRA_LIBRARIES} gfortran gomp)
   #target_link_libraries(${argon} PRIVATE nblib libdftbplus libdftd3 libmudpack gfortran gomp)

    set_target_properties(${methane} PROPERTIES LINK_FLAGS
      "-Wl,--whole-archive ${GMX_QMMM_DFTBPLUS_LIB}/libdftbplus.a -Wl,--no-whole-archive")
    target_link_libraries(${methane} PRIVATE nblib libdftbplus ${LINEAR_ALGEBRA_LIBRARIES} gfortran gomp)
else()
    target_link_libraries(${argon} PRIVATE nblib)
    target_link_libraries(${methane} PRIVATE nblib)
endif()

if(BUILD_TESTING)
    add_test(NAME NbLibSamplesTestArgon COMMAND ${argon})
    add_test(NAME NbLibSamplesTestMethaneWater COMMAND ${methane})
    add_dependencies(nblib-tests ${argon} ${methane})
endif()
