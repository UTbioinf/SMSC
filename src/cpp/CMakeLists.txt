cmake_minimum_required(VERSION 3.6)

project(smsc_proj)

find_package(LEMON REQUIRED)
if(LEMON_FOUND)
    include_directories(${LEMON_INCLUDE_DIRS})
endif()

message("-- [ROOT]: Install path [${CMAKE_INSTALL_PREFIX}]")
#set(install_submodule ON CACHE BOOL "Instlal libs" FORCE)
add_subdirectory(mummer-lib)
add_subdirectory(common-bio)
add_subdirectory(smsc)
