include_directories ("${PROJECT_BINARY_DIR}/"
"${PROJECT_SOURCE_DIR}/src/pfs" "${CMAKE_CURRENT_SOURCE_DIR}")

link_directories("${PROJECT_SOURCE_DIR}/src/pfs")


add_executable(linearhdr pfshdrraw.cpp responses.cpp linearhdr.cpp nrutil.cpp "${GETOPT_OBJECT}")
target_link_libraries(linearhdr pfs)
# install (TARGETS ${TRG2} DESTINATION bin)



