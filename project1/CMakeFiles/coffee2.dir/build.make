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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ekakshkataria/git/smooth/project1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ekakshkataria/git/smooth/project1

# Include any dependencies generated for this target.
include CMakeFiles/coffee2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/coffee2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/coffee2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/coffee2.dir/flags.make

CMakeFiles/coffee2.dir/coffee2.cpp.o: CMakeFiles/coffee2.dir/flags.make
CMakeFiles/coffee2.dir/coffee2.cpp.o: coffee2.cpp
CMakeFiles/coffee2.dir/coffee2.cpp.o: CMakeFiles/coffee2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ekakshkataria/git/smooth/project1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/coffee2.dir/coffee2.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/coffee2.dir/coffee2.cpp.o -MF CMakeFiles/coffee2.dir/coffee2.cpp.o.d -o CMakeFiles/coffee2.dir/coffee2.cpp.o -c /Users/ekakshkataria/git/smooth/project1/coffee2.cpp

CMakeFiles/coffee2.dir/coffee2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/coffee2.dir/coffee2.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ekakshkataria/git/smooth/project1/coffee2.cpp > CMakeFiles/coffee2.dir/coffee2.cpp.i

CMakeFiles/coffee2.dir/coffee2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/coffee2.dir/coffee2.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ekakshkataria/git/smooth/project1/coffee2.cpp -o CMakeFiles/coffee2.dir/coffee2.cpp.s

# Object files for target coffee2
coffee2_OBJECTS = \
"CMakeFiles/coffee2.dir/coffee2.cpp.o"

# External object files for target coffee2
coffee2_EXTERNAL_OBJECTS =

coffee2: CMakeFiles/coffee2.dir/coffee2.cpp.o
coffee2: CMakeFiles/coffee2.dir/build.make
coffee2: /Users/ekakshkataria/git/smooth/software/lib/libsmooth.a
coffee2: /Users/ekakshkataria/git/commonutils/software/lib/libmsu_commonutils.a
coffee2: /opt/homebrew/lib/libgsl.dylib
coffee2: /opt/homebrew/lib/libgslcblas.dylib
coffee2: CMakeFiles/coffee2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ekakshkataria/git/smooth/project1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable coffee2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/coffee2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/coffee2.dir/build: coffee2
.PHONY : CMakeFiles/coffee2.dir/build

CMakeFiles/coffee2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/coffee2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/coffee2.dir/clean

CMakeFiles/coffee2.dir/depend:
	cd /Users/ekakshkataria/git/smooth/project1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ekakshkataria/git/smooth/project1 /Users/ekakshkataria/git/smooth/project1 /Users/ekakshkataria/git/smooth/project1 /Users/ekakshkataria/git/smooth/project1 /Users/ekakshkataria/git/smooth/project1/CMakeFiles/coffee2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/coffee2.dir/depend

