# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/quirren/git/smooth/bin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/quirren/git/smooth/bin

# Include any dependencies generated for this target.
include CMakeFiles/simplex.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/simplex.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/simplex.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simplex.dir/flags.make

CMakeFiles/simplex.dir/codes/simplex.cc.o: CMakeFiles/simplex.dir/flags.make
CMakeFiles/simplex.dir/codes/simplex.cc.o: codes/simplex.cc
CMakeFiles/simplex.dir/codes/simplex.cc.o: CMakeFiles/simplex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/simplex.dir/codes/simplex.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/simplex.dir/codes/simplex.cc.o -MF CMakeFiles/simplex.dir/codes/simplex.cc.o.d -o CMakeFiles/simplex.dir/codes/simplex.cc.o -c /home/quirren/git/smooth/bin/codes/simplex.cc

CMakeFiles/simplex.dir/codes/simplex.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simplex.dir/codes/simplex.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/bin/codes/simplex.cc > CMakeFiles/simplex.dir/codes/simplex.cc.i

CMakeFiles/simplex.dir/codes/simplex.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simplex.dir/codes/simplex.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/bin/codes/simplex.cc -o CMakeFiles/simplex.dir/codes/simplex.cc.s

# Object files for target simplex
simplex_OBJECTS = \
"CMakeFiles/simplex.dir/codes/simplex.cc.o"

# External object files for target simplex
simplex_EXTERNAL_OBJECTS =

simplex: CMakeFiles/simplex.dir/codes/simplex.cc.o
simplex: CMakeFiles/simplex.dir/build.make
simplex: /home/quirren/git/smooth/software/lib/libsmooth.a
simplex: /home/quirren/git/commonutils/software/lib/libmsu_commonutils.a
simplex: /usr/lib/x86_64-linux-gnu/libgsl.so
simplex: /usr/lib/x86_64-linux-gnu/libgslcblas.so
simplex: CMakeFiles/simplex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/quirren/git/smooth/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable simplex"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simplex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simplex.dir/build: simplex
.PHONY : CMakeFiles/simplex.dir/build

CMakeFiles/simplex.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simplex.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simplex.dir/clean

CMakeFiles/simplex.dir/depend:
	cd /home/quirren/git/smooth/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quirren/git/smooth/bin /home/quirren/git/smooth/bin /home/quirren/git/smooth/bin /home/quirren/git/smooth/bin /home/quirren/git/smooth/bin/CMakeFiles/simplex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/simplex.dir/depend

