#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2011- The GROMACS Authors
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

add_executable(template template.cpp)
if (MSVC)
    gmx_target_warning_suppression(template /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION_DOUBLE_TO_REAL)
    gmx_target_warning_suppression(template /wd4267 HAS_NO_MSVC_LOSSY_CONVERSION_SIZE_T_TO_INT)
endif()
# This should be removable once object libraries can directly use target_link_libraries
# with CMake 3.12, #3290
target_include_directories(template SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
target_link_libraries(template PRIVATE libgromacs legacy_modules legacy_api ${GMX_EXE_LINKER_FLAGS})

# Ensure that libsycl is properly linked when libgromacs is a static
# library.
#
# TODO The second predicate is used in release-2025 branch to
# maximize stability of the default shared-library build
# configuration. After merging to main branch, remove the second
# predicate, because it may be useful and seems unlikely to cause
# harm.
if(GMX_GPU_SYCL AND NOT BUILD_SHARED_LIBS)
    # There are no source files to gmx that need to compile as
    # SYCL, so none are provided here.
    add_sycl_to_target(TARGET template)
endif()

if (GMX_QMMM_DFTBPLUS)
    LINK_DIRECTORIES(GMX_QMMM_DFTBPLUS_LIB)
    add_library(libdftbplus STATIC IMPORTED)
    set_target_properties(libdftbplus PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libdftbplus.a)
   #add_library(libxmlf90 STATIC IMPORTED)
   #set_target_properties(libxmlf90 PROPERTIES IMPORTED_LOCATION ${GMX_QMMM_DFTBPLUS_LIB}/libxmlf90.a)
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

#   add_library(gfortran-7.2 SHARED IMPORTED)
#   set_target_properties(gfortran-7.2 PROPERTIES IMPORTED_LOCATION /usr/local/src/gcc-7.2.0-release/lib64/libgfortran.so.4)
    target_link_libraries(template libdftbplus gfortran gomp)
   #target_link_libraries(template libdftbplus libdftd3 libmudpack gfortran gomp)
   #target_link_libraries(template libdftbplus libxmlf90 libdftd3 gfortran gomp)
endif()

set(DOCUMENTATION_HTTP_URL_BASE
    https://manual.gromacs.org/nightly/index.html)
if (SOURCE_IS_SOURCE_DISTRIBUTION)
    set(DOCUMENTATION_HTTP_URL_BASE
        http://manual.gromacs.org/documentation/${GMX_VERSION_STRING})
endif()
configure_file(README.cmakein README @ONLY)

install(FILES CMakeLists.txt.template
        DESTINATION ${GMX_INSTALL_GMXDATADIR}/template
        RENAME CMakeLists.txt
        COMPONENT development)

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/README template.cpp Makefile.pkg
        DESTINATION ${GMX_INSTALL_GMXDATADIR}/template
        COMPONENT development)

install(FILES cmake/FindGROMACS.cmake
        DESTINATION ${GMX_INSTALL_GMXDATADIR}/template/cmake
        COMPONENT development)
