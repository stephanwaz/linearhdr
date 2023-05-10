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

configure_file(src/linearhdr_make_list.py ${CMAKE_INSTALL_PREFIX}/bin/linearhdr_make_list COPYONLY FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
configure_file(src/linearhdr_calibrate.py ${CMAKE_INSTALL_PREFIX}/bin/linearhdr_calibrate COPYONLY FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
configure_file(README.rst ${CMAKE_INSTALL_PREFIX}/bin/linearhdr_README.rst COPYONLY)
