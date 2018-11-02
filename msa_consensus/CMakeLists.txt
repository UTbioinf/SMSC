cmake_minimum_required(VERSION 3.6)


project (msa_consensus_lib)

option(install_submodule "Install common_bio" OFF)

include_directories(${msa_lib_SOURCE_DIR}/..)
add_library(msa_consensus msa_consensus.cpp)
target_link_libraries(msa_consensus msa_lib parasail)

if(install_submodule)
    install(TARGETS msa_consensus DESTINATION lib)
    install(FILES msa_consensus.h DESTINATION include/msa_consensus)
endif()
