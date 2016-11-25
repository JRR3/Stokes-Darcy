# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/javier/dealii-8.3.0/examples/darcy_p

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/javier/dealii-8.3.0/examples/darcy_p

# Utility rule file for strip_comments.

# Include the progress variables for this target.
include CMakeFiles/strip_comments.dir/progress.make

CMakeFiles/strip_comments:
	$(CMAKE_COMMAND) -E cmake_progress_report /home/javier/dealii-8.3.0/examples/darcy_p/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "strip comments"
	/usr/bin/perl -pi -e 's#^[ \t]*//.*\n##g;' convergence_rates.cc flux_function_2.cc stokes_exact_solution_2.cc interfacial_pressure_function_2.cc darcy_exact_solution_2.cc darcy_K_inverse.cc plot_i_data.cc map_linker.cc parallel_map_linker.cc mesh_data.cc build_lambda_flux_p_maps.cc my_l2_norm.cc sdd_simple_parallel.cc

strip_comments: CMakeFiles/strip_comments
strip_comments: CMakeFiles/strip_comments.dir/build.make
.PHONY : strip_comments

# Rule to build all files generated by this target.
CMakeFiles/strip_comments.dir/build: strip_comments
.PHONY : CMakeFiles/strip_comments.dir/build

CMakeFiles/strip_comments.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/strip_comments.dir/cmake_clean.cmake
.PHONY : CMakeFiles/strip_comments.dir/clean

CMakeFiles/strip_comments.dir/depend:
	cd /home/javier/dealii-8.3.0/examples/darcy_p && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/javier/dealii-8.3.0/examples/darcy_p /home/javier/dealii-8.3.0/examples/darcy_p /home/javier/dealii-8.3.0/examples/darcy_p /home/javier/dealii-8.3.0/examples/darcy_p /home/javier/dealii-8.3.0/examples/darcy_p/CMakeFiles/strip_comments.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/strip_comments.dir/depend

