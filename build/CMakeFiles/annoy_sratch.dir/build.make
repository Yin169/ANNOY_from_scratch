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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.4/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/yincheangng/worksapce/C++/annoy_from_sratch

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/yincheangng/worksapce/C++/annoy_from_sratch/build

# Include any dependencies generated for this target.
include CMakeFiles/annoy_sratch.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/annoy_sratch.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/annoy_sratch.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/annoy_sratch.dir/flags.make

CMakeFiles/annoy_sratch.dir/main.cpp.o: CMakeFiles/annoy_sratch.dir/flags.make
CMakeFiles/annoy_sratch.dir/main.cpp.o: /Users/yincheangng/worksapce/C++/annoy_from_sratch/main.cpp
CMakeFiles/annoy_sratch.dir/main.cpp.o: CMakeFiles/annoy_sratch.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yincheangng/worksapce/C++/annoy_from_sratch/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/annoy_sratch.dir/main.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/annoy_sratch.dir/main.cpp.o -MF CMakeFiles/annoy_sratch.dir/main.cpp.o.d -o CMakeFiles/annoy_sratch.dir/main.cpp.o -c /Users/yincheangng/worksapce/C++/annoy_from_sratch/main.cpp

CMakeFiles/annoy_sratch.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/annoy_sratch.dir/main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yincheangng/worksapce/C++/annoy_from_sratch/main.cpp > CMakeFiles/annoy_sratch.dir/main.cpp.i

CMakeFiles/annoy_sratch.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/annoy_sratch.dir/main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yincheangng/worksapce/C++/annoy_from_sratch/main.cpp -o CMakeFiles/annoy_sratch.dir/main.cpp.s

# Object files for target annoy_sratch
annoy_sratch_OBJECTS = \
"CMakeFiles/annoy_sratch.dir/main.cpp.o"

# External object files for target annoy_sratch
annoy_sratch_EXTERNAL_OBJECTS =

annoy_sratch: CMakeFiles/annoy_sratch.dir/main.cpp.o
annoy_sratch: CMakeFiles/annoy_sratch.dir/build.make
annoy_sratch: CMakeFiles/annoy_sratch.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/yincheangng/worksapce/C++/annoy_from_sratch/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable annoy_sratch"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/annoy_sratch.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/annoy_sratch.dir/build: annoy_sratch
.PHONY : CMakeFiles/annoy_sratch.dir/build

CMakeFiles/annoy_sratch.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/annoy_sratch.dir/cmake_clean.cmake
.PHONY : CMakeFiles/annoy_sratch.dir/clean

CMakeFiles/annoy_sratch.dir/depend:
	cd /Users/yincheangng/worksapce/C++/annoy_from_sratch/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/yincheangng/worksapce/C++/annoy_from_sratch /Users/yincheangng/worksapce/C++/annoy_from_sratch /Users/yincheangng/worksapce/C++/annoy_from_sratch/build /Users/yincheangng/worksapce/C++/annoy_from_sratch/build /Users/yincheangng/worksapce/C++/annoy_from_sratch/build/CMakeFiles/annoy_sratch.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/annoy_sratch.dir/depend

