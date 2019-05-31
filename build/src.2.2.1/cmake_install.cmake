# Install script for directory: /afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/src.2.2.1

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/build/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/blasrm5_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/blasrm5_filter")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/blasrm5_filter"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/build/src.2.2.1/blasrm5_filter")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/blasrm5_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/blasrm5_filter")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/afs/crc.nd.edu/user/s/szhu3/loonlocal_intel/linuxbrew/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/blasrm5_filter")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/smsc" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/smsc")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/smsc"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/build/src.2.2.1/smsc")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/smsc" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/smsc")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/afs/crc.nd.edu/user/s/szhu3/loonlocal_intel/linuxbrew/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/smsc")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/src.2.2.1/run_nucmer.sh"
    "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/src.2.2.1/config.sh"
    "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/src.2.2.1/run_mhap.sh"
    "/afs/crc.nd.edu/user/s/szhu3/Private/github/UTbioinf/SMSC/src.2.2.1/run_muscle.sh"
    )
endif()

