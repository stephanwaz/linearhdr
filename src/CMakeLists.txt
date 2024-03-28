add_subdirectory (pfs)


include_directories ("${PROJECT_BINARY_DIR}/"
        "${PROJECT_SOURCE_DIR}/src/pfs" "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/fileformat")

link_directories("${PROJECT_SOURCE_DIR}/src/pfs")


include_directories(${TIFF_INCLUDE_DIR})

add_executable(linearhdr linearhdr_main.cpp linearhdr.cpp dht_demosaic.cpp fileformat/rgbeio.cpp fileformat/hdrtiffio.cpp "${GETOPT_OBJECT}")
target_link_libraries(linearhdr pfs ${TIFF_LIBRARY})

add_executable(convertcolor convertcolor.cpp fileformat/rgbeio.cpp linearhdr.cpp dht_demosaic.cpp "${GETOPT_OBJECT}")
target_link_libraries(convertcolor pfs)

message(STATUS ${CMAKE_INSTALL_PREFIX})

install (TARGETS linearhdr DESTINATION bin)
install (TARGETS convertcolor DESTINATION bin)
