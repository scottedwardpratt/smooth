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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.22.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.22.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/scottpratt/git/smooth/scottrun

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/scottpratt/git/smooth/scottrun

# Utility rule file for extern_commonutils.

# Include any custom commands dependencies for this target.
include CMakeFiles/extern_commonutils.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/extern_commonutils.dir/progress.make

CMakeFiles/extern_commonutils:
	cd /Users/scottpratt/git/commonutils/software && make

extern_commonutils: CMakeFiles/extern_commonutils
extern_commonutils: CMakeFiles/extern_commonutils.dir/build.make
.PHONY : extern_commonutils

# Rule to build all files generated by this target.
CMakeFiles/extern_commonutils.dir/build: extern_commonutils
.PHONY : CMakeFiles/extern_commonutils.dir/build

CMakeFiles/extern_commonutils.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/extern_commonutils.dir/cmake_clean.cmake
.PHONY : CMakeFiles/extern_commonutils.dir/clean

CMakeFiles/extern_commonutils.dir/depend:
	cd /Users/scottpratt/git/smooth/scottrun && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun/CMakeFiles/extern_commonutils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/extern_commonutils.dir/depend

