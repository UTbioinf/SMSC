# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/loon/.linuxbrew/Cellar/cmake/3.11.2/bin/cmake

# The command to remove a file.
RM = /home/loon/.linuxbrew/Cellar/cmake/3.11.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon

# Include any dependencies generated for this target.
include test/CMakeFiles/bpgraph_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/bpgraph_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/bpgraph_test.dir/flags.make

test/CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.o: test/CMakeFiles/bpgraph_test.dir/flags.make
test/CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.o: /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1/test/bpgraph_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.o"
	cd /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.o -c /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1/test/bpgraph_test.cc

test/CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.i"
	cd /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1/test/bpgraph_test.cc > CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.i

test/CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.s"
	cd /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1/test/bpgraph_test.cc -o CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.s

# Object files for target bpgraph_test
bpgraph_test_OBJECTS = \
"CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.o"

# External object files for target bpgraph_test
bpgraph_test_EXTERNAL_OBJECTS =

test/bpgraph_test: test/CMakeFiles/bpgraph_test.dir/bpgraph_test.cc.o
test/bpgraph_test: test/CMakeFiles/bpgraph_test.dir/build.make
test/bpgraph_test: lemon/libemon.a
test/bpgraph_test: test/CMakeFiles/bpgraph_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bpgraph_test"
	cd /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bpgraph_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/bpgraph_test.dir/build: test/bpgraph_test

.PHONY : test/CMakeFiles/bpgraph_test.dir/build

test/CMakeFiles/bpgraph_test.dir/clean:
	cd /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test && $(CMAKE_COMMAND) -P CMakeFiles/bpgraph_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/bpgraph_test.dir/clean

test/CMakeFiles/bpgraph_test.dir/depend:
	cd /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1 /home/loon/Loon/github/UTbioinf/SMSC/depends/lemon-1.3.1/test /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test /home/loon/Loon/github/UTbioinf/SMSC/depends/build_lemon/test/CMakeFiles/bpgraph_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/bpgraph_test.dir/depend

