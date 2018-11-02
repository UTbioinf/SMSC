cmake_minimum_required(VERSION 3.6)

project(smsc_proj VERSION 2.2.1)

message("-- [ROOT]: Install path [${CMAKE_INSTALL_PREFIX}]")
add_subdirectory(mummer-lib)
add_subdirectory(common-bio)
#add_subdirectory(msa_lib) # replaced by muscle
#add_subdirectory(msa_consensus) # put into src/, but will be integrated with a modified muscle and resume this
add_subdirectory(src.${PROJECT_VERSION})
