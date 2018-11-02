cmake_minimum_required(VERSION 3.6)

project(common_bio)

option(install_submodule "Install common_bio" OFF)

include_directories(${mummer_api_SOURCE_DIR} ${mummer_libbase_SOURCE_DIR} ${mummer_streesrc_SOURCE_DIR})
add_library(common-bio Aligner.cpp Consensus.cpp Multiseq.cpp Progress.cpp baseheader.cpp nucleotide.cpp util.cpp)
target_link_libraries(common-bio mummer base stree)

if(install_submodule)
    set(LIBCOMMON_BIO_HEADER Aligner.h Consensus.h Multiseq.h Progress.h baseheader.h multi-array.h nucleotide.h util.h)
    install(TARGETS common-bio DESTINATION lib)
    install(FILES ${LIBCOMMON_BIO_HEADER} DESTINATION include/common-bio)
endif()
