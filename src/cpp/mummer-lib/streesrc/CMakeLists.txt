cmake_minimum_required(VERSION 3.6)

project(mummer_streesrc)

option(install_submodule "Install libstree.a" OFF)

include_directories(${mummer_libbase_SOURCE_DIR})
add_library(stree access.c depthtab.c iterator.c overmax.c stree.c addleafcount.c dfs.c linkloc.c oversucc.c construct.c ex2leav.c loc.c scanpref.c)
target_link_libraries(stree base)
target_compile_definitions(stree PUBLIC STREEHUGE)

if(install_submodule)
    set(LIBSTREE_HEADERS streeacc.h streehuge.h streemac.h streesmall.h symboldef.h streedef.h streelarge.h streeproto.h streetyp.h)
    install(TARGETS stree DESTINATION lib)
    install(FILES ${LIBSTREE_HEADERS} DESTINATION include/mummer)
endif()
