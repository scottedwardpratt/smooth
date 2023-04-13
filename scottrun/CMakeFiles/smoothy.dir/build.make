# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/scottpratt/git/smooth/scottrun

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/scottpratt/git/smooth/scottrun

# Include any dependencies generated for this target.
include CMakeFiles/smoothy.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/smoothy.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/smoothy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/smoothy.dir/flags.make

CMakeFiles/smoothy.dir/smoothmain.cc.o: CMakeFiles/smoothy.dir/flags.make
CMakeFiles/smoothy.dir/smoothmain.cc.o: smoothmain.cc
CMakeFiles/smoothy.dir/smoothmain.cc.o: CMakeFiles/smoothy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/scottpratt/git/smooth/scottrun/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/smoothy.dir/smoothmain.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smoothy.dir/smoothmain.cc.o -MF CMakeFiles/smoothy.dir/smoothmain.cc.o.d -o CMakeFiles/smoothy.dir/smoothmain.cc.o -c /Users/scottpratt/git/smooth/scottrun/smoothmain.cc

CMakeFiles/smoothy.dir/smoothmain.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smoothy.dir/smoothmain.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/scottpratt/git/smooth/scottrun/smoothmain.cc > CMakeFiles/smoothy.dir/smoothmain.cc.i

CMakeFiles/smoothy.dir/smoothmain.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smoothy.dir/smoothmain.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/scottpratt/git/smooth/scottrun/smoothmain.cc -o CMakeFiles/smoothy.dir/smoothmain.cc.s

# Object files for target smoothy
smoothy_OBJECTS = \
"CMakeFiles/smoothy.dir/smoothmain.cc.o"

# External object files for target smoothy
smoothy_EXTERNAL_OBJECTS =

smoothy: CMakeFiles/smoothy.dir/smoothmain.cc.o
smoothy: CMakeFiles/smoothy.dir/build.make
smoothy: /Users/scottpratt/git/smooth/software/lib/libsmooth.a
smoothy: /Users/scottpratt/git/coral/software/lib/libmsu_coral.a
smoothy: /Users/scottpratt/git/commonutils/software/lib/libmsu_commonutils.a
smoothy: /usr/local/Cellar/gsl/2.7.1/lib/libgsl.dylib
smoothy: /usr/local/Cellar/gsl/2.7.1/lib/libgslcblas.dylib
smoothy: CMakeFiles/smoothy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/scottpratt/git/smooth/scottrun/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable smoothy"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/smoothy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/smoothy.dir/build: smoothy
.PHONY : CMakeFiles/smoothy.dir/build

CMakeFiles/smoothy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/smoothy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/smoothy.dir/clean

CMakeFiles/smoothy.dir/depend:
	cd /Users/scottpratt/git/smooth/scottrun && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun/CMakeFiles/smoothy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/smoothy.dir/depend

