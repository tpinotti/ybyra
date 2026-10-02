# puera - visualizing y-chromosome samples on trees
# Copyright (C) 2025 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@sund.ku.dk>
# University of Copenhagen, Globe Institute, Section for GeoGenetics
# Oster Voldgade 5-7, 1350 Copenhagen K, Denmark

# --------------------------------------------------------------------------------------------------
#   CMake Init
# --------------------------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

# Fun with colors!
if (NOT WIN32)
    string(ASCII 27 Esc)

    set(ColorBlue   "${Esc}[34m") # Build information.
    set(ColorGreen  "${Esc}[32m") # Optional additions.
    set(ColorYellow "${Esc}[33m") # Tipps for the user.
    set(ColorRed    "${Esc}[31m") # Warnings and Errors.
    set(ColorEnd    "${Esc}[0m")
endif()

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set (CMAKE_DISABLE_SOURCE_CHANGES  ON)

# Need to make sure that AppleClang is actually recognized and not mistakingly identifyed
# as normal Clang. See https://cmake.org/cmake/help/v3.0/policy/CMP0025.html
# This is important for the test cmake script, where we need to set some warnings for clang,
# but not for apple clang.
if( NOT ( CMAKE_VERSION VERSION_LESS 3.0 ))
    # Appleclang compatibility
    cmake_policy(SET CMP0025 NEW)
    # Support for :: aliasing, used in new way to include OpenMP
    cmake_policy(SET CMP0028 NEW)
endif()

# --------------------------------------------------------------------------------------------------
#   Project Definition and Requirements
# --------------------------------------------------------------------------------------------------

project ( puera CXX )

# --------------------------------------------------------------------------------------------------
#   Build Options
# --------------------------------------------------------------------------------------------------

# The build type is determined by the standard cmake variable CMAKE_BUILD_TYPE.
# Thus, if you want to change the build type, call `cmake .. -DCMAKE_BUILD_TYPE=DEBUG`
IF(CMAKE_BUILD_TYPE STREQUAL "")
    set (CMAKE_BUILD_TYPE RELEASE)
ENDIF()

# Print build type. Useful for debugging user issues.
message (STATUS "${ColorBlue}puera build type: ${CMAKE_BUILD_TYPE}${ColorEnd}")

# Store everything in `bin`.
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )

# Option to enable building a fully static version, that is, statical linking of system libraries.
# This works with GCC, but might cause trouble with OpenMP and/or Clang, so it is disabled by
# default. Most users who compile on their own probaby won't need the fully portable version.
# Genesis on the other hand is always linked statically, in order to be able to move around the
# final binary at least on the same system.
option ( PUERA_BUILD_STATIC  "Build the static version of the executable." OFF )

if( PUERA_BUILD_STATIC AND ( NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" ))
    string(ASCII 27 Esc)
    message( STATUS "${Esc}[31m"
        "Compiler is not GCC. This might cause trouble with statically linking system libraries "
        "on Linux systems. If the compilation fails with some linking error, try a different compiler, "
        "or try to not build statically by deactivating the CMake option 'PUERA_BUILD_STATIC'. "
    "${Esc}[0m" )
endif()

# --------------------------------------------------------------------------------------------------
#   Compiler and Linker Options
# --------------------------------------------------------------------------------------------------

if (MSVC OR (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC"))
    # MinGW/MSVC do not support much. Also, the above string comparison in the condition does
    # not work there, hence the extra `MSVC` check...
    set (WARN_FLAGS "-Wall")
else ()
    # For other operating systems, set the warn flags to a very high level - except unknown pragmas,
    # as this causes needless warnings with OpenMP and other pragma-based techniques.
    set (WARN_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -pedantic -pedantic-errors")

    # Furthermore, in debug mode, we also want all warnings to be errors, so that we do not miss them
    # In release, this might hinder users to successfully build genesis if any warnings remain
    # or occur only on their compiler, so let's avoid trouble for them.
    IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
        set (WARN_FLAGS "${WARN_FLAGS} -Werror")
    ENDIF()
endif ()

set (CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++17 ${WARN_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG   "-O2 -DDEBUG -g -ggdb3 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" )

IF(PUERA_BUILD_STATIC)
    # Use different ways of telling CMake to build static, just to be sure.
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    link_libraries("-static")
ENDIF()
message( STATUS "Static linking of system libraries: ${PUERA_BUILD_STATIC}")

# Genesis uses LTO, and so do we.
if( NOT ( CMAKE_VERSION VERSION_LESS 3.9 ))
    cmake_policy(SET CMP0069 NEW)
endif()

# --------------------------------------------------------------------------------------------------
#   Dependencies
# --------------------------------------------------------------------------------------------------

# Genesis commit to download if the submodule is not checked out (e.g., when downloaded as a zip).
# Set by tools/cmake/update_dependencies.sh to the checked out submodule commit; do not edit manually.
SET( genesis_COMMIT_HASH "b4e4f4101791a4f3abde76d64b0afe4424e43a79" ) #genesis_COMMIT_HASH#

# Use the genesis submodule if present, otherwise download the pinned commit from GitHub.
# We only use FetchContent for the download (the non-existing SOURCE_SUBDIR prevents it from adding
# the project), and add it ourselves, so that the variables genesis exports reach our scope.
message (STATUS "Looking for genesis")
set( GENESIS_SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/libs/genesis )
if( EXISTS "${GENESIS_SUBMODULE_DIR}/CMakeLists.txt" )
    message( STATUS "Using genesis submodule" )
    set( GENESIS_SOURCE_DIR "${GENESIS_SUBMODULE_DIR}" )
else()
    message( STATUS "${ColorBlue}Using genesis from GitHub @ ${genesis_COMMIT_HASH}${ColorEnd}" )
    include( FetchContent )
    FetchContent_Declare(
        genesis
        URL "https://github.com/lczech/genesis/archive/${genesis_COMMIT_HASH}.zip"
        SOURCE_SUBDIR "no-add-subdirectory"
    )
    FetchContent_MakeAvailable( genesis )
    set( GENESIS_SOURCE_DIR "${genesis_SOURCE_DIR}" )
endif()

# Policy CMP0042: Enable MACOSX_RPATH
if (APPLE)
   set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)

# Force to set the options used for genesis. Bit cumbersome in cmake...
# We use C++17, so we deactivate the genesis auto detection of the highest available C++ standard.
# We do not need OpenMP or htslib, so deactivate them to keep the build simple.
set(GENESIS_CPP_STANDARD "17" CACHE STRING "" FORCE)
set(GENESIS_USE_OPENMP OFF CACHE BOOL "" FORCE)
set(GENESIS_USE_HTSLIB OFF CACHE BOOL "" FORCE)

# Add Genesis as dependency.
add_subdirectory( ${GENESIS_SOURCE_DIR} ${CMAKE_BINARY_DIR}/genesis )

# Use everything that Genesis exports, just to be sure that we use the same setup.
add_definitions( ${GENESIS_DEFINITIONS} )
include_directories( ${GENESIS_INCLUDE_DIR} )
set( CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${GENESIS_C_FLAGS}")
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )

# ------------------------------------------------------------------------------
#   Sources
# ------------------------------------------------------------------------------

# Find all source files. The CMake documentation does not recommend this,
# but for our purposes, this is good enough for now.
file (GLOB_RECURSE puera_sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Add our own code. The main is only one cpp file, the rest are headers included from there.
include_directories( "src" )
add_executable( puera ${puera_sources} )

# Link it against Genesis and all its dependencies, and against the CLI11 target provided by genesis.
target_link_libraries ( puera ${GENESIS_LINK_LIBRARIES} CLI11 )

# We have to use LTO for puera as well for some compilers to properly link, I think.
if(GENESIS_HAS_LTO)
    message(STATUS "Building with LTO/IPO support")
    set_property(TARGET puera PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "Building without LTO/IPO support")
endif()

message( STATUS "CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}" )
message( STATUS "GENESIS_LINK_LIBRARIES ${GENESIS_LINK_LIBRARIES}" )
