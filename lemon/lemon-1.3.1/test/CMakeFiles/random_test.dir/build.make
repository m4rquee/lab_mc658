# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.12.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.12.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1

# Include any dependencies generated for this target.
include test/CMakeFiles/random_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/random_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/random_test.dir/flags.make

test/CMakeFiles/random_test.dir/random_test.cc.o: test/CMakeFiles/random_test.dir/flags.make
test/CMakeFiles/random_test.dir/random_test.cc.o: test/random_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/random_test.dir/random_test.cc.o"
	cd /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/random_test.dir/random_test.cc.o -c /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test/random_test.cc

test/CMakeFiles/random_test.dir/random_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/random_test.dir/random_test.cc.i"
	cd /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test/random_test.cc > CMakeFiles/random_test.dir/random_test.cc.i

test/CMakeFiles/random_test.dir/random_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/random_test.dir/random_test.cc.s"
	cd /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test/random_test.cc -o CMakeFiles/random_test.dir/random_test.cc.s

# Object files for target random_test
random_test_OBJECTS = \
"CMakeFiles/random_test.dir/random_test.cc.o"

# External object files for target random_test
random_test_EXTERNAL_OBJECTS =

test/random_test: test/CMakeFiles/random_test.dir/random_test.cc.o
test/random_test: test/CMakeFiles/random_test.dir/build.make
test/random_test: lemon/libemon.a
test/random_test: /usr/local/lib/libglpk.dylib
test/random_test: test/CMakeFiles/random_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable random_test"
	cd /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/random_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/random_test.dir/build: test/random_test

.PHONY : test/CMakeFiles/random_test.dir/build

test/CMakeFiles/random_test.dir/clean:
	cd /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test && $(CMAKE_COMMAND) -P CMakeFiles/random_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/random_test.dir/clean

test/CMakeFiles/random_test.dir/depend:
	cd /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1 /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1 /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test /Users/fkm/Dropbox/tiziu/programming/branch_and_cut_and_connectivity/connectivity/lab_mc658/lemon/lemon-1.3.1/test/CMakeFiles/random_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/random_test.dir/depend

