cmake_minimum_required(VERSION 3.6)

project(mummer_libbasedir)

option(install_submodule "Install libbase.a" OFF)

add_library(base cleanMUMcand.c multiseq.c safescpy.c space.c clock.c mapfile.c procopt.c seterror.c)

if(install_submodule)
    set(LIBBASE_HEADERS args.h debugdef.h intbits.h multidef.h protodef.h visible.h arraydef.h errordef.h megabytes.h mumcand.h spacedef.h chardef.h fopen.h minmax.h optdesc.h types.h)
    install(TARGETS base DESTINATION lib)
    install(FILES ${LIBBASE_HEADERS} DESTINATION include/mummer)
endif()
