
include_directories ("${PROJECT_BINARY_DIR}/")
if( NOT HAS_GETOPT )
	include_directories ("${GETOPT_INCLUDE}")
endif( NOT HAS_GETOPT )
add_library(pfs ${LIB_TYPE} colorspace.cpp pfs.cpp pfsutils.cpp array2d.h pfs.h "${GETOPT_OBJECT}")

# SOVERSION changes only when the library API changes, so it may be
# different from the package version
set_target_properties(pfs PROPERTIES VERSION 2.0.0 SOVERSION 2)



