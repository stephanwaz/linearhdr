cmake_minimum_required (VERSION 2.8.8)

set(CMAKE_INSTALL_PREFIX /Users/stephenwasilewski/code/linearhdr/build)

project(linearhdr)

include_directories(/opt/local/include/)

set(VERSION_MAJOR 2)
set(VERSION_MINOR 1)

set(PKG_DATADIR "${CMAKE_INSTALL_PREFIX}/share/pfstools")
configure_file (
        "${PROJECT_SOURCE_DIR}/cmake_config.h.in"
        "${PROJECT_BINARY_DIR}/config.h"
)

add_subdirectory (src) 
