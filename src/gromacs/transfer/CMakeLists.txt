#
# Simulate electron and exciton transfer.
#
# Dept. Theoretical Chemical Biology
# Karlsruhe Institute of Technology
#
# Copyright (c) 2020 onwards

# Set up the module library
add_library(transfer INTERFACE)
file(GLOB TRANSFER_SOURCES *.cpp)
list(REMOVE_ITEM TRANSFER_SOURCES dftb-obsolete.cpp leftovers.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${TRANSFER_SOURCES} PARENT_SCOPE)

set(TRANSFER_SOURCES ${TRANSFER_SOURCES} PARENT_SCOPE)

## Source files have the following private module dependencies.
#target_link_libraries(transfer PRIVATE
##                      gmxlib
##                      math
##                      mdtypes
##                      tng_io
#                      )

# Public interface for modules, including dependencies and interfaces
target_include_directories(transfer INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(transfer INTERFACE
                      legacy_api
                      )

# TODO: when transfer is an OBJECT target
#target_link_libraries(transfer PUBLIC legacy_api)
#target_link_libraries(transfer PRIVATE common)

# Module dependencies
# mdlib interfaces convey transitive dependence on these modules.
#target_link_libraries(transfer PUBLIC
target_link_libraries(transfer INTERFACE
                      utility
                      )
# Source files have the following private module dependencies.
#target_link_libraries(transfer PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(transfer PRIVATE legacy_modules)

