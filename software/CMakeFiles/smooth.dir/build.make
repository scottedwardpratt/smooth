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
CMAKE_SOURCE_DIR = /home/quirren/git/smooth/software

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/quirren/git/smooth/software

# Include any dependencies generated for this target.
include CMakeFiles/smooth.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/smooth.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/smooth.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/smooth.dir/flags.make

CMakeFiles/smooth.dir/src/emulator.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/emulator.cc.o: src/emulator.cc
CMakeFiles/smooth.dir/src/emulator.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/smooth.dir/src/emulator.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/emulator.cc.o -MF CMakeFiles/smooth.dir/src/emulator.cc.o.d -o CMakeFiles/smooth.dir/src/emulator.cc.o -c /home/quirren/git/smooth/software/src/emulator.cc

CMakeFiles/smooth.dir/src/emulator.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/emulator.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/emulator.cc > CMakeFiles/smooth.dir/src/emulator.cc.i

CMakeFiles/smooth.dir/src/emulator.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/emulator.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/emulator.cc -o CMakeFiles/smooth.dir/src/emulator.cc.s

CMakeFiles/smooth.dir/src/latincube.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/latincube.cc.o: src/latincube.cc
CMakeFiles/smooth.dir/src/latincube.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/smooth.dir/src/latincube.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/latincube.cc.o -MF CMakeFiles/smooth.dir/src/latincube.cc.o.d -o CMakeFiles/smooth.dir/src/latincube.cc.o -c /home/quirren/git/smooth/software/src/latincube.cc

CMakeFiles/smooth.dir/src/latincube.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/latincube.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/latincube.cc > CMakeFiles/smooth.dir/src/latincube.cc.i

CMakeFiles/smooth.dir/src/latincube.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/latincube.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/latincube.cc -o CMakeFiles/smooth.dir/src/latincube.cc.s

CMakeFiles/smooth.dir/src/master.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/master.cc.o: src/master.cc
CMakeFiles/smooth.dir/src/master.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/smooth.dir/src/master.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/master.cc.o -MF CMakeFiles/smooth.dir/src/master.cc.o.d -o CMakeFiles/smooth.dir/src/master.cc.o -c /home/quirren/git/smooth/software/src/master.cc

CMakeFiles/smooth.dir/src/master.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/master.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/master.cc > CMakeFiles/smooth.dir/src/master.cc.i

CMakeFiles/smooth.dir/src/master.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/master.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/master.cc -o CMakeFiles/smooth.dir/src/master.cc.s

CMakeFiles/smooth.dir/src/modelparinfo.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/modelparinfo.cc.o: src/modelparinfo.cc
CMakeFiles/smooth.dir/src/modelparinfo.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/smooth.dir/src/modelparinfo.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/modelparinfo.cc.o -MF CMakeFiles/smooth.dir/src/modelparinfo.cc.o.d -o CMakeFiles/smooth.dir/src/modelparinfo.cc.o -c /home/quirren/git/smooth/software/src/modelparinfo.cc

CMakeFiles/smooth.dir/src/modelparinfo.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/modelparinfo.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/modelparinfo.cc > CMakeFiles/smooth.dir/src/modelparinfo.cc.i

CMakeFiles/smooth.dir/src/modelparinfo.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/modelparinfo.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/modelparinfo.cc -o CMakeFiles/smooth.dir/src/modelparinfo.cc.s

CMakeFiles/smooth.dir/src/observableinfo.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/observableinfo.cc.o: src/observableinfo.cc
CMakeFiles/smooth.dir/src/observableinfo.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/smooth.dir/src/observableinfo.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/observableinfo.cc.o -MF CMakeFiles/smooth.dir/src/observableinfo.cc.o.d -o CMakeFiles/smooth.dir/src/observableinfo.cc.o -c /home/quirren/git/smooth/software/src/observableinfo.cc

CMakeFiles/smooth.dir/src/observableinfo.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/observableinfo.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/observableinfo.cc > CMakeFiles/smooth.dir/src/observableinfo.cc.i

CMakeFiles/smooth.dir/src/observableinfo.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/observableinfo.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/observableinfo.cc -o CMakeFiles/smooth.dir/src/observableinfo.cc.s

CMakeFiles/smooth.dir/src/real.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/real.cc.o: src/real.cc
CMakeFiles/smooth.dir/src/real.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/smooth.dir/src/real.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/real.cc.o -MF CMakeFiles/smooth.dir/src/real.cc.o.d -o CMakeFiles/smooth.dir/src/real.cc.o -c /home/quirren/git/smooth/software/src/real.cc

CMakeFiles/smooth.dir/src/real.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/real.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/real.cc > CMakeFiles/smooth.dir/src/real.cc.i

CMakeFiles/smooth.dir/src/real.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/real.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/real.cc -o CMakeFiles/smooth.dir/src/real.cc.s

CMakeFiles/smooth.dir/src/scorecard.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/scorecard.cc.o: src/scorecard.cc
CMakeFiles/smooth.dir/src/scorecard.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/smooth.dir/src/scorecard.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/scorecard.cc.o -MF CMakeFiles/smooth.dir/src/scorecard.cc.o.d -o CMakeFiles/smooth.dir/src/scorecard.cc.o -c /home/quirren/git/smooth/software/src/scorecard.cc

CMakeFiles/smooth.dir/src/scorecard.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/scorecard.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/scorecard.cc > CMakeFiles/smooth.dir/src/scorecard.cc.i

CMakeFiles/smooth.dir/src/scorecard.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/scorecard.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/scorecard.cc -o CMakeFiles/smooth.dir/src/scorecard.cc.s

CMakeFiles/smooth.dir/src/simplex.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/simplex.cc.o: src/simplex.cc
CMakeFiles/smooth.dir/src/simplex.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/smooth.dir/src/simplex.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/simplex.cc.o -MF CMakeFiles/smooth.dir/src/simplex.cc.o.d -o CMakeFiles/smooth.dir/src/simplex.cc.o -c /home/quirren/git/smooth/software/src/simplex.cc

CMakeFiles/smooth.dir/src/simplex.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/simplex.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/simplex.cc > CMakeFiles/smooth.dir/src/simplex.cc.i

CMakeFiles/smooth.dir/src/simplex.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/simplex.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/simplex.cc -o CMakeFiles/smooth.dir/src/simplex.cc.s

CMakeFiles/smooth.dir/src/smooth.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/smooth.cc.o: src/smooth.cc
CMakeFiles/smooth.dir/src/smooth.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/smooth.dir/src/smooth.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/smooth.cc.o -MF CMakeFiles/smooth.dir/src/smooth.cc.o.d -o CMakeFiles/smooth.dir/src/smooth.cc.o -c /home/quirren/git/smooth/software/src/smooth.cc

CMakeFiles/smooth.dir/src/smooth.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/smooth.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/smooth.cc > CMakeFiles/smooth.dir/src/smooth.cc.i

CMakeFiles/smooth.dir/src/smooth.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/smooth.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/smooth.cc -o CMakeFiles/smooth.dir/src/smooth.cc.s

CMakeFiles/smooth.dir/src/traininginfo.cc.o: CMakeFiles/smooth.dir/flags.make
CMakeFiles/smooth.dir/src/traininginfo.cc.o: src/traininginfo.cc
CMakeFiles/smooth.dir/src/traininginfo.cc.o: CMakeFiles/smooth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/smooth.dir/src/traininginfo.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/smooth.dir/src/traininginfo.cc.o -MF CMakeFiles/smooth.dir/src/traininginfo.cc.o.d -o CMakeFiles/smooth.dir/src/traininginfo.cc.o -c /home/quirren/git/smooth/software/src/traininginfo.cc

CMakeFiles/smooth.dir/src/traininginfo.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/smooth.dir/src/traininginfo.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quirren/git/smooth/software/src/traininginfo.cc > CMakeFiles/smooth.dir/src/traininginfo.cc.i

CMakeFiles/smooth.dir/src/traininginfo.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/smooth.dir/src/traininginfo.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quirren/git/smooth/software/src/traininginfo.cc -o CMakeFiles/smooth.dir/src/traininginfo.cc.s

# Object files for target smooth
smooth_OBJECTS = \
"CMakeFiles/smooth.dir/src/emulator.cc.o" \
"CMakeFiles/smooth.dir/src/latincube.cc.o" \
"CMakeFiles/smooth.dir/src/master.cc.o" \
"CMakeFiles/smooth.dir/src/modelparinfo.cc.o" \
"CMakeFiles/smooth.dir/src/observableinfo.cc.o" \
"CMakeFiles/smooth.dir/src/real.cc.o" \
"CMakeFiles/smooth.dir/src/scorecard.cc.o" \
"CMakeFiles/smooth.dir/src/simplex.cc.o" \
"CMakeFiles/smooth.dir/src/smooth.cc.o" \
"CMakeFiles/smooth.dir/src/traininginfo.cc.o"

# External object files for target smooth
smooth_EXTERNAL_OBJECTS =

lib/libsmooth.a: CMakeFiles/smooth.dir/src/emulator.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/latincube.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/master.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/modelparinfo.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/observableinfo.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/real.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/scorecard.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/simplex.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/smooth.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/src/traininginfo.cc.o
lib/libsmooth.a: CMakeFiles/smooth.dir/build.make
lib/libsmooth.a: CMakeFiles/smooth.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/quirren/git/smooth/software/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX static library lib/libsmooth.a"
	$(CMAKE_COMMAND) -P CMakeFiles/smooth.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/smooth.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/smooth.dir/build: lib/libsmooth.a
.PHONY : CMakeFiles/smooth.dir/build

CMakeFiles/smooth.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/smooth.dir/cmake_clean.cmake
.PHONY : CMakeFiles/smooth.dir/clean

CMakeFiles/smooth.dir/depend:
	cd /home/quirren/git/smooth/software && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quirren/git/smooth/software /home/quirren/git/smooth/software /home/quirren/git/smooth/software /home/quirren/git/smooth/software /home/quirren/git/smooth/software/CMakeFiles/smooth.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/smooth.dir/depend

