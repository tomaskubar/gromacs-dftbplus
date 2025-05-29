#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014 by the GROMACS development team.
# Copyright (c) 2015,2016,2018,2019,2020, by the GROMACS development team, led by
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

file(GLOB MDRUN_SOURCES mdrun/*.cpp)
# make an "object library" that we can re-use for multiple targets
add_library(mdrun_objlib OBJECT ${MDRUN_SOURCES})
gmx_target_compile_options(mdrun_objlib)
target_compile_definitions(mdrun_objlib PRIVATE HAVE_CONFIG_H)
target_include_directories(mdrun_objlib SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
# Should be possible to remove this when resolving #3290
target_include_directories(mdrun_objlib SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
target_link_libraries(mdrun_objlib PRIVATE common)
target_link_libraries(mdrun_objlib PRIVATE legacy_api)
# TODO: Explicitly link specific modules.
target_link_libraries(mdrun_objlib PRIVATE legacy_modules)

if(GMX_FAHCORE)
    # The lack of a real source file here alongside the object library
    # may break some generators, according to CMake documentation. If
    # so, we can consider adding some dummy file to make it work.
    add_library(fahcore $<TARGET_OBJECTS:mdrun_objlib>)
    target_link_libraries(fahcore PRIVATE ${GMX_COMMON_LIBRARIES} legacy_api)
elseif(GMX_BUILD_MDRUN_ONLY)
    add_executable(mdrun-only $<TARGET_OBJECTS:mdrun_objlib> mdrun_main.cpp)
    gmx_target_compile_options(mdrun-only)
    target_include_directories(mdrun-only SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
    target_compile_definitions(mdrun-only PRIVATE HAVE_CONFIG_H)
    target_link_libraries(mdrun-only PRIVATE
                          common
                          legacy_modules
                          libgromacs
                          ${GMX_COMMON_LIBRARIES}
                          ${GMX_EXE_LINKER_FLAGS}
                          )
    set(BINARY_NAME "mdrun${GMX_BINARY_SUFFIX}")
    set_target_properties(mdrun-only PROPERTIES
        OUTPUT_NAME "${BINARY_NAME}")
    install(TARGETS mdrun-only DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT mdrun-only)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gmx-completion-${BINARY_NAME}.bash
         "complete -o nospace -F _gmx_mdrun_compl ${BINARY_NAME}")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gmx-completion-${BINARY_NAME}.bash
            DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
else()
    file(GLOB GMX_MAIN_SOURCES gmx.cpp legacymodules.cpp)
    if(GMX_X11)
        file(GLOB VIEW_SOURCES view/*.cpp)
    else()
        file(GLOB VIEW_SOURCES view/view.cpp)
    endif()
    add_library(view_objlib OBJECT ${VIEW_SOURCES})
    gmx_target_compile_options(view_objlib)
    target_compile_definitions(view_objlib PRIVATE HAVE_CONFIG_H)
    target_include_directories(view_objlib SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
    target_link_libraries(view_objlib PRIVATE common legacy_api)
    # TODO: Explicitly link specific modules.
    target_link_libraries(view_objlib PRIVATE legacy_modules)
    add_library(gmx_objlib OBJECT ${GMX_MAIN_SOURCES})
    target_link_libraries(gmx_objlib PRIVATE common legacy_api)
    # TODO: Explicitly link specific modules.
    target_link_libraries(gmx_objlib PRIVATE legacy_modules)
    target_include_directories(gmx_objlib SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
    target_include_directories(gmx_objlib SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)
    add_executable(gmx
        $<TARGET_OBJECTS:gmx_objlib>
        $<TARGET_OBJECTS:mdrun_objlib>
        $<TARGET_OBJECTS:view_objlib>)
    add_executable(Gromacs::gmx ALIAS gmx)
    gmx_target_compile_options(gmx)
    target_compile_definitions(gmx PRIVATE HAVE_CONFIG_H)
    target_link_libraries(gmx PRIVATE
                          common
                          libgromacs
                          ${GMX_COMMON_LIBRARIES}
                          ${GMX_EXE_LINKER_FLAGS})
    if(GMX_X11)
        target_link_libraries(gmx ${X11_LIBRARIES})
    endif()
    set(BINARY_NAME "gmx${GMX_BINARY_SUFFIX}")
    set_target_properties(gmx PROPERTIES
        OUTPUT_NAME "${BINARY_NAME}")
    # Add the executable to the INSTALL target and to the libgromacs EXPORT configuration.
    # The GROMACS installation then provides an importable (e.g. find_package(GROMACS))
    # Gromacs::gmx target, which may be used as the executable COMMAND argument for functions like
    # add_custom_command() or add_test(). The full path to the installed binary may retrieved
    # from the LOCATION property or with a $<TARGET_FILE:Gromacs:gmx> generator expression
    # (where generator expressions are allowed).
    install(TARGETS gmx
            EXPORT libgromacs
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    ########################
    # Completion generation

    include(gmxCustomCommandUtilities)

    set(COMPLETION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/completion)
    # Using GMX_BUILD_HELP here is somewhat confusing, but the conditions when
    # this can be done are exactly the same (ability to run the compiled
    # binaries).
    if (GMX_BUILD_HELP)
        gmx_add_custom_output_target(completion OUTPUT STAMP
            COMMAND ${CMAKE_COMMAND}
                -D GMX_EXECUTABLE=$<TARGET_FILE:gmx>
                -D ERRORS_ARE_FATAL=${GMX_BUILD_HELP_FORCE}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/BuildCompletions.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            DEPENDS gmx ${CMAKE_CURRENT_SOURCE_DIR}/BuildCompletions.cmake
            COMMENT "Generating command-line completions for programs")
        set_target_properties(completion PROPERTIES EXCLUDE_FROM_ALL OFF)
        set_directory_properties(PROPERTIES
            ADDITIONAL_MAKE_CLEAN_FILES "completion")
        set(COMPLETION_DIR ${CMAKE_CURRENT_BINARY_DIR}/completion)
    endif()
    if (SOURCE_IS_SOURCE_DISTRIBUTION OR GMX_BUILD_HELP)
        install(DIRECTORY ${COMPLETION_DIR}/
                DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime OPTIONAL)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/gmx-completion-${BINARY_NAME}.bash
             "complete -o nospace -F _gmx_compl ${BINARY_NAME}")
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gmx-completion-${BINARY_NAME}.bash
                DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
    endif()
    gmx_cpack_add_generated_source_directory(completion)

    if(BUILD_TESTING)
        add_subdirectory(mdrun/tests)
    endif()
endif()

if (GMX_QMMM_DFTBPLUS OR GMX_QMMM_DFTBPLUS_PYTORCH OR GMX_QMMM_DFTBPLUS_TENSORFLOW)
    LINK_DIRECTORIES(GMX_QMMM_DFTBPLUS_LIB)
    
    # Find OpenBLAS explicitly
    find_library(OPENBLAS_LIBRARY openblas PATHS /usr/local/run/OpenBLAS-0.3.20-be/lib)
    if(NOT OPENBLAS_LIBRARY)
        message(FATAL_ERROR "OpenBLAS required by DFTB+ was not found")
    endif()
    message(STATUS "Found OpenBLAS: ${OPENBLAS_LIBRARY}")
    
    # Find ARPACK library
    find_library(ARPACK_LIBRARY arpack PATHS /usr/local/run/arpack-ng-3.7.0-openblas/lib)
    if(NOT ARPACK_LIBRARY)
        message(FATAL_ERROR "ARPACK required by DFTB+ was not found")
    endif()
    message(STATUS "Found ARPACK: ${ARPACK_LIBRARY}")
    
    # Create a unified dependency list with ARPACK included
    set(DFTBPLUS_DEPS libdftbplus libdftd3 libmudpack gfortran gomp ${OPENBLAS_LIBRARY} ${ARPACK_LIBRARY})

    # Add dependencies to both libgromacs and executables
    target_link_libraries(libgromacs PRIVATE ${DFTBPLUS_DEPS})

    # Rest of your existing DFTB+ setup...
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

    if(GMX_BUILD_MDRUN_ONLY)
        target_link_libraries(mdrun PRIVATE ${DFTBPLUS_DEPS})
       #target_link_libraries(mdrun libdftbplus libxmlf90 libdftd3 gfortran gomp)
    else()
        target_link_libraries(gmx PRIVATE ${DFTBPLUS_DEPS})
       #target_link_libraries(gmx libdftbplus libxmlf90 libdftd3 gfortran gomp)
    endif()
endif()

if("${GMX_QMMM_PROGRAM}" STREQUAL "TENSORFLOW" OR "${GMX_QMMM_PROGRAM}" STREQUAL "DFTBPLUS_TENSORFLOW")
    find_library(TENSORFLOW_LIB
        NAMES libtensorflow.so
        PATHS /usr/local/lib $ENV{TENSORFLOW_ROOT}
        HINTS $ENV{TENSORFLOW_ROOT}/lib ${CMAKE_PREFIX_PATH}
        DOC "Find tensorflow library libtensorflow.so"
        REQUIRED
    )
    # find_library(TENSORFLOW_FRAMEWORK_LIB
    #     NAMES libtensorflow_framework.so
    #     PATHS /usr/local/lib $ENV{TENSORFLOW_ROOT}/lib
    #     HINTS $ENV{TENSORFLOW_ROOT} ${CMAKE_PREFIX_PATH}
    #     DOC "Find tensorflow framework library libtensorflow_framework.so"
    #     REQUIRED
    # )
    find_library(JSONCPP_LIB
        NAMES libjsoncpp.so libjsoncpp.a 
        PATHS /usr/local/lib $ENV{JSONCPP_ROOT}/lib
        HINTS $ENV{JSONCPP_ROOT} ${CMAKE_PREFIX_PATH}
        DOC "Find jsoncpp library libjsoncpp.so"
        REQUIRED
    )

    if(TENSORFLOW_LIB)
        message(STATUS "Looking for tensorflow at ${TENSORFLOW_LIB}")
    else()
        message(WARNING "TensorFlow library not found. Please set the TENSORFLOW_ROOT environment variable or specify the library path manually using -DTENSORFLOW_LIB=<path_to_libtensorflow.so>.")
        message(FATAL_ERROR "tensorflow library not found.")
    endif()

    if (TENSORFLOW_FRAMEWORK_LIB)
        message(STATUS "Looking for tensorflow framework library at ${TENSORFLOW_FRAMEWORK_LIB}")
    else()
        message(WARNING "TensorFlow framework library not found. Please set the TENSORFLOW_ROOT environment variable or specify the library path manually using -DTENSORFLOW_FRAMEWORK_LIB=<path_to_libtensorflow_framework.so>.")
        message(FATAL_ERROR "tensorflow framework library not found.")
    endif()
    if (JSONCPP_LIB)
        message(STATUS "Looking for jsoncpp library at ${JSONCPP_LIB}")
    else()
        message(WARNING "jsoncpp library not found. Please set the JSONCPP_ROOT environment variable or specify the library path manually using -DJSONCPP_LIB=<path_to_libjsoncpp.so>.")
        message(FATAL_ERROR "jsoncpp library not found.")
    endif()

    add_library(libtensorflow STATIC IMPORTED)
    set_target_properties(libtensorflow PROPERTIES IMPORTED_LOCATION ${TENSORFLOW_LIB})
    add_library(libtensorflow_framework STATIC IMPORTED)
    set_target_properties(libtensorflow_framework PROPERTIES IMPORTED_LOCATION ${TENSORFLOW_FRAMEWORK_LIB})
    # add_library(libjsoncpp STATIC IMPORTED)
    # set_target_properties(libjsoncpp PROPERTIES IMPORTED_LOCATION ${JSONCPP_LIB})

    # Link TensorFlow to both libgromacs and gmx
    target_link_libraries(libgromacs PRIVATE "${TENSORFLOW_LIB}" "${JSONCPP_LIB}")
    #target_include_directories(libgromacs PRIVATE "${TENSORFLOW_INCLUDE_DIR}")
    target_link_libraries(gmx PRIVATE "${TENSORFLOW_LIB}" "${JSONCPP_LIB}")
    #target_include_directories(gmx PRIVATE "${TENSORFLOW_INCLUDE_DIR}")
endif()

if ("${GMX_QMMM_PROGRAM}" STREQUAL "PYTORCH" OR "${GMX_QMMM_PROGRAM}" STREQUAL "DFTBPLUS_PYTORCH") 
    find_package(Torch REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

    add_library(libtorch STATIC IMPORTED)
    set_target_properties(libtorch PROPERTIES IMPORTED_LOCATION "${TORCH_LIBRARIES}")
    
    target_link_libraries(gmx PRIVATE "${TORCH_LIBRARIES}")
    target_include_directories(gmx PRIVATE "${TORCH_INCLUDE_DIRS}")
    
    set_property(TARGET libtorch PROPERTY CXX_STANDARD 17)
    message(STATUS "Torch found (include: ${TORCH_INCLUDE_DIRS}, lib: ${TORCH_LIBRARIES})")
endif()
