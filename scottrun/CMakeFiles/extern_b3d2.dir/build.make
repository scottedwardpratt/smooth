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

# Utility rule file for extern_b3d2.

# Include any custom commands dependencies for this target.
include CMakeFiles/extern_b3d2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/extern_b3d2.dir/progress.make

CMakeFiles/extern_b3d2:
	cd /Users/scottpratt/git/hp/alice/software/b3d2 && make

extern_b3d2: CMakeFiles/extern_b3d2
extern_b3d2: CMakeFiles/extern_b3d2.dir/build.make
.PHONY : extern_b3d2

# Rule to build all files generated by this target.
CMakeFiles/extern_b3d2.dir/build: extern_b3d2
.PHONY : CMakeFiles/extern_b3d2.dir/build

CMakeFiles/extern_b3d2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/extern_b3d2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/extern_b3d2.dir/clean

CMakeFiles/extern_b3d2.dir/depend:
	cd /Users/scottpratt/git/smooth/scottrun && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun /Users/scottpratt/git/smooth/scottrun/CMakeFiles/extern_b3d2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/extern_b3d2.dir/depend

