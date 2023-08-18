#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
#
# \author Victor Holanda <victor.holanda@cscs.ch>
# \author Joe Jordan <ejjordan@kth.se>
# \author Prashanth Kanduri <kanduri@cscs.ch>
# \author Sebastian Keller <keller@cscs.ch>
#

add_executable(argon-forces-integration "")
add_library(libtensorflow SHARED IMPORTED)
set_target_properties(libtensorflow PROPERTIES IMPORTED_LOCATION /usr/local/lib/libtensorflow.so)

target_sources(argon-forces-integration
        PRIVATE
        argon-forces-integration.cpp
        )

target_link_libraries(argon-forces-integration PRIVATE nblib)
install(
        TARGETS
        argon-forces-integration
        DESTINATION
        ${CMAKE_INSTALL_BINDIR}
)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

if (GMX_QMMM_DFTBPLUS)
    LINK_DIRECTORIES(GMX_QMMM_DFTBPLUS_LIB)
    add_library(libdftbplus STATIC IMPORTED)
    set_target_properties(libdftbplus PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libdftbplus.a)
    add_library(libdftd3 STATIC IMPORTED)
    set_target_properties(libdftd3 PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libdftd3.a)
    add_library(libmudpack STATIC IMPORTED)
    set_target_properties(libmudpack PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libmudpack.a)

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

    target_link_libraries(argon-forces-integration PRIVATE libdftbplus libdftd3 libmudpack gfortran gomp)
endif()

target_link_libraries(libtensorflow)