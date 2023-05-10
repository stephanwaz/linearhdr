include_directories ("${PROJECT_BINARY_DIR}/"
"${PROJECT_SOURCE_DIR}/src/pfs" "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/fileformat")

link_directories("${PROJECT_SOURCE_DIR}/src/pfs")

include_directories(${NETPBM_INCLUDE_DIR})


if( NETPBM_FOUND )
    add_executable(linearhdr linearhdr_main.cpp responses.cpp linearhdr.cpp nrutil.cpp ../fileformat/ppmio.cpp ../fileformat/rgbeio.cpp "${GETOPT_OBJECT}")
    target_link_libraries(linearhdr pfs ${NETPBM_LIBRARIES})
else()
    add_executable(linearhdr linearhdr_main.cpp responses.cpp linearhdr.cpp nrutil.cpp ../fileformat/rgbeio.cpp "${GETOPT_OBJECT}")
    target_link_libraries(linearhdr pfs)
endif()

message(STATUS ${CMAKE_INSTALL_PREFIX})

install (TARGETS linearhdr DESTINATION bin)