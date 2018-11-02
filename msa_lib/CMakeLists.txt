cmake_minimum_required(VERSION 3.6)


project ("msa_lib")

option(install_submodule "Install common_bio" OFF)

add_library(msa_lib msa_bias.cpp msa_ecalc.cpp msa_lib.cpp)
target_link_libraries(msa_lib parasail)

if(install_submodule)
    install(TARGETS msa_lib DESTINATION lib)
    install(FILES msa_base.h msa_bias.h msa_ecalc.h msa_lib.h DESTINATION include/msa_lib)
endif()
