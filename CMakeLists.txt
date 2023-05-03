cmake_minimum_required (VERSION 2.8.8)


project(linearhdr)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/build)

include_directories(/opt/local/include/)

set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

OPTION(WITH_NetPBM "Use NetPBM library" ON)

if( WITH_NetPBM )
#    set(CMAKE_FIND_DEBUG_MODE TRUE)
    find_package(NETPBM)
#    set(CMAKE_FIND_DEBUG_MODE FALSE)
    if( NOT NETPBM_FOUND )
        MESSAGE( STATUS "NetPBM not found." )
    endif( NOT NETPBM_FOUND )

else( WITH_NetPBM )
    set( NETPBM_FOUND OFF )
endif( WITH_NetPBM )

set(PKG_DATADIR "${CMAKE_INSTALL_PREFIX}/share/pfstools")
configure_file (
        "${PROJECT_SOURCE_DIR}/cmake_config.h.in"
        "${PROJECT_BINARY_DIR}/config.h"
)

add_subdirectory (src) 
