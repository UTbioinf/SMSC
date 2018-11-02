cmake_minimum_required(VERSION 3.6)

project(mummer_api)

option(install_submodule "Install mummer_api" OFF)

add_subdirectory(libbasedir)
add_subdirectory(streesrc)

include_directories(${mummer_libbase_SOURCE_DIR} ${mummer_streesrc_SOURCE_DIR})
add_definitions(-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE)
add_library(mummer delta-filter.cpp mgaps.cpp postnuc.cpp sw_align.cpp delta.cpp nucmer.cpp prenuc.cpp tigrinc.cpp findmaxmat.c)
target_link_libraries(mummer base stree)

if(install_submodule)
    set(LIBMUMMER_HEADER delta-filter.h maxmatdef.h nucmer.h prenuc.h sw_alignscore.h delta.h mgaps.h postnuc.h sw_align.h tigrinc.h)
    install(TARGETS mummer DESTINATION lib)
    install(FILES ${LIBMUMMER_HEADER} DESTINATION include/mummer)
endif()
