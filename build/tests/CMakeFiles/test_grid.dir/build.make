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
CMAKE_COMMAND = "C:/Program Files/CMake/bin/cmake.exe"

# The command to remove a file.
RM = "C:/Program Files/CMake/bin/cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:/cpp/Sci/Flow

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:/cpp/Sci/Flow/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/test_grid.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/test_grid.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/test_grid.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/test_grid.dir/flags.make

tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj: tests/CMakeFiles/test_grid.dir/flags.make
tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj: tests/CMakeFiles/test_grid.dir/includes_CXX.rsp
tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj: ../tests/testGridFromFile.cpp
tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj: tests/CMakeFiles/test_grid.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:/cpp/Sci/Flow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj"
	cd D:/cpp/Sci/Flow/build/tests && C:/Strawberry/c/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj -MF CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj.d -o CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj -c D:/cpp/Sci/Flow/tests/testGridFromFile.cpp

tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_grid.dir/testGridFromFile.cpp.i"
	cd D:/cpp/Sci/Flow/build/tests && C:/Strawberry/c/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:/cpp/Sci/Flow/tests/testGridFromFile.cpp > CMakeFiles/test_grid.dir/testGridFromFile.cpp.i

tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_grid.dir/testGridFromFile.cpp.s"
	cd D:/cpp/Sci/Flow/build/tests && C:/Strawberry/c/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:/cpp/Sci/Flow/tests/testGridFromFile.cpp -o CMakeFiles/test_grid.dir/testGridFromFile.cpp.s

# Object files for target test_grid
test_grid_OBJECTS = \
"CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj"

# External object files for target test_grid
test_grid_EXTERNAL_OBJECTS =

tests/test_grid.exe: tests/CMakeFiles/test_grid.dir/testGridFromFile.cpp.obj
tests/test_grid.exe: tests/CMakeFiles/test_grid.dir/build.make
tests/test_grid.exe: src/libvovka_library.a
tests/test_grid.exe: _deps/fmtlib-build/libfmtd.a
tests/test_grid.exe: tests/CMakeFiles/test_grid.dir/linklibs.rsp
tests/test_grid.exe: tests/CMakeFiles/test_grid.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:/cpp/Sci/Flow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_grid.exe"
	cd D:/cpp/Sci/Flow/build/tests && "C:/Program Files/CMake/bin/cmake.exe" -E rm -f CMakeFiles/test_grid.dir/objects.a
	cd D:/cpp/Sci/Flow/build/tests && C:/Strawberry/c/bin/ar.exe qc CMakeFiles/test_grid.dir/objects.a @CMakeFiles/test_grid.dir/objects1.rsp
	cd D:/cpp/Sci/Flow/build/tests && C:/Strawberry/c/bin/g++.exe -g -Wl,--whole-archive CMakeFiles/test_grid.dir/objects.a -Wl,--no-whole-archive -o test_grid.exe -Wl,--out-implib,libtest_grid.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/test_grid.dir/linklibs.rsp

# Rule to build all files generated by this target.
tests/CMakeFiles/test_grid.dir/build: tests/test_grid.exe
.PHONY : tests/CMakeFiles/test_grid.dir/build

tests/CMakeFiles/test_grid.dir/clean:
	cd D:/cpp/Sci/Flow/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_grid.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/test_grid.dir/clean

tests/CMakeFiles/test_grid.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" D:/cpp/Sci/Flow D:/cpp/Sci/Flow/tests D:/cpp/Sci/Flow/build D:/cpp/Sci/Flow/build/tests D:/cpp/Sci/Flow/build/tests/CMakeFiles/test_grid.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/test_grid.dir/depend

