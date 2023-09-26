cmake_minimum_required (VERSION 3.0.0)


project(linearhdr)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/build CACHE PATH "installation destination (parent of bin)" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)

include_directories(/opt/local/include/)
include( CheckCXXSourceCompiles )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# ======== Check for getopt, which is missing on Windows ========

check_cxx_source_compiles( "#include \"getopt.h\" \n int main() { return 0; }" HAS_GETOPT )
if( NOT HAS_GETOPT )
    MESSAGE( STATUS "getopt not found, a replacement needs to be compiled. " )
    set( GETOPT_INCLUDE "${PROJECT_SOURCE_DIR}/src/getopt" )
    set( GETOPT_OBJECT "$<TARGET_OBJECTS:getopt>" )
else( NOT HAS_GETOPT )
    MESSAGE( STATUS "getopt found. " )
    set( GETOPT_INCLUDE ".") # work around, cannot pass empty string
    set( GETOPT_OBJECT )
endif( NOT HAS_GETOPT )


find_package(TIFF)

set(PKG_DATADIR "${CMAKE_INSTALL_PREFIX}/share/pfstools")
configure_file (
        "${PROJECT_SOURCE_DIR}/cmake_config.h.in"
        "${PROJECT_BINARY_DIR}/config.h"
)


set(LIBRAW_PATH ${CMAKE_CURRENT_SOURCE_DIR}/LibRaw)
file(TO_CMAKE_PATH "${LIBRAW_PATH}" LIBRAW_PATH)

# ==================================================================================================
# Library version info extraction

file(READ ${LIBRAW_PATH}/libraw/libraw_version.h _libraw_version_content)

# API version strings
string(REGEX MATCH "#define LIBRAW_MAJOR_VERSION[ \t]*([0-9]*)\n" _version_major_match ${_libraw_version_content})
set(RAW_LIB_MAJOR_VERSION "${CMAKE_MATCH_1}")

string(REGEX MATCH "#define LIBRAW_MINOR_VERSION[ \t]*([0-9]*)\n" _version_minor_match ${_libraw_version_content})
set(RAW_LIB_MINOR_VERSION "${CMAKE_MATCH_1}")

string(REGEX MATCH "#define LIBRAW_PATCH_VERSION[ \t]*([0-9]*)\n" _version_patch_match ${_libraw_version_content})
set(RAW_LIB_PATCH_VERSION "${CMAKE_MATCH_1}")

# Set env. variables accordinly.
set(RAW_LIB_VERSION_STRING    "${RAW_LIB_MAJOR_VERSION}.${RAW_LIB_MINOR_VERSION}.${RAW_LIB_PATCH_VERSION}")

message(STATUS "LibRaw string version: ${RAW_LIB_VERSION_STRING}")


if(WIN32 AND NOT DEFINED CMAKE_DEBUG_POSTFIX)
    set(CMAKE_DEBUG_POSTFIX "d")
endif()

# To prevent warnings from M$ compiler
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_ATL_SECURE_NO_WARNINGS)
    add_definitions(-D_AFX_SECURE_NO_WARNINGS)
endif()


# Math library check

if(NOT WIN32)
    FIND_LIBRARY(MATH_LIBRARY m)
endif()

# Put the include dirs which are in the source or build tree
# before all other include dirs, so the headers in the sources
# are preferred over the already installed ones
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/
        ${LIBRAW_PATH}/
        )


file(GLOB_RECURSE libraw_LIB_SRCS CONFIGURE_DEPENDS "${LIBRAW_PATH}/src/*.cpp")

# Exclude placeholder (stub) implementations
file(GLOB_RECURSE exclude_libraw_LIB_SRCS CONFIGURE_DEPENDS "${LIBRAW_PATH}/src/*_ph.cpp")
list(REMOVE_ITEM libraw_LIB_SRCS ${exclude_libraw_LIB_SRCS})


# linearhdr patch to dcraw_process
list(TRANSFORM libraw_LIB_SRCS REPLACE LibRaw/src/postprocessing/dcraw_process.cpp src/dcraw_process.cpp)

add_library(raw ${libraw_LIB_SRCS})
add_library(libraw::libraw ALIAS raw)
target_compile_definitions(raw PRIVATE LIBRAW_NOTHREADS)


add_subdirectory(src)

configure_file(src/bin/extract_data.sh ${CMAKE_INSTALL_PREFIX}/bin/linearhdr_extract COPYONLY FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# add a small macro so that this is a bit cleaner
macro(LIBRAW_BUILD_SAMPLES)
    set(_filename ${ARGV0})
    set(_rawlib ${ARGV1})
    set(_basepath ${ARGV2})
    string(REPLACE "." ";" _temp ${_filename})
    list(GET _temp 0 _target)

    set(${_target}_SRCS ${_basepath}/${_filename})

    add_executable(${_target} ${${_target}_SRCS})
    target_compile_options(${_target} PRIVATE -w)

    target_link_libraries(${_target} PRIVATE ${_rawlib})

    if(${_rawlib} MATCHES "raw_r")
        target_link_libraries(${_target} PUBLIC ${PTHREADS_LIBRARY})
    endif()

    install (TARGETS ${_target} DESTINATION bin)

endmacro(LIBRAW_BUILD_SAMPLES)

LIBRAW_BUILD_SAMPLES(rawconvert.cpp raw src)


