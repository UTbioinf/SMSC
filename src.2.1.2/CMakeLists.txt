cmake_minimum_required(VERSION 3.6)

project(smsc)

include_directories(${CMAKE_SOURCE_DIR} ${mummer_api_SOURCE_DIR} ${mummer_libbase_SOURCE_DIR} ${mummer_streesrc_SOURCE_DIR})

add_executable(blasrm5_filter blasrm5_filter.cpp)

add_executable(smsc smsc.cpp mummer_api.cpp threadpool.c path_cover.cpp)
set_source_files_properties(threadpool.c PROPERTIES COMPILE_FLAGS -O3 -fPIC -D_REENTRANT -Wall -pedantic)
target_link_libraries(smsc common-bio mummer base stree emon pthread)

install(TARGETS blasrm5_filter smsc DESTINATION bin)
install(PROGRAMS run_nucmer.sh DESTINATION bin)
