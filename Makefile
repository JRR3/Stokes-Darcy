# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/javier/dealii-8.3.0/examples/darcy_p/CMakeFiles /home/javier/dealii-8.3.0/examples/darcy_p/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/javier/dealii-8.3.0/examples/darcy_p/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named debug

# Build rule for target.
debug: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 debug
.PHONY : debug

# fast build rule for target.
debug/fast:
	$(MAKE) -f CMakeFiles/debug.dir/build.make CMakeFiles/debug.dir/build
.PHONY : debug/fast

#=============================================================================
# Target rules for targets named distclean

# Build rule for target.
distclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 distclean
.PHONY : distclean

# fast build rule for target.
distclean/fast:
	$(MAKE) -f CMakeFiles/distclean.dir/build.make CMakeFiles/distclean.dir/build
.PHONY : distclean/fast

#=============================================================================
# Target rules for targets named info

# Build rule for target.
info: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 info
.PHONY : info

# fast build rule for target.
info/fast:
	$(MAKE) -f CMakeFiles/info.dir/build.make CMakeFiles/info.dir/build
.PHONY : info/fast

#=============================================================================
# Target rules for targets named release

# Build rule for target.
release: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 release
.PHONY : release

# fast build rule for target.
release/fast:
	$(MAKE) -f CMakeFiles/release.dir/build.make CMakeFiles/release.dir/build
.PHONY : release/fast

#=============================================================================
# Target rules for targets named run

# Build rule for target.
run: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 run
.PHONY : run

# fast build rule for target.
run/fast:
	$(MAKE) -f CMakeFiles/run.dir/build.make CMakeFiles/run.dir/build
.PHONY : run/fast

#=============================================================================
# Target rules for targets named runclean

# Build rule for target.
runclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 runclean
.PHONY : runclean

# fast build rule for target.
runclean/fast:
	$(MAKE) -f CMakeFiles/runclean.dir/build.make CMakeFiles/runclean.dir/build
.PHONY : runclean/fast

#=============================================================================
# Target rules for targets named sdd_simple_parallel

# Build rule for target.
sdd_simple_parallel: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sdd_simple_parallel
.PHONY : sdd_simple_parallel

# fast build rule for target.
sdd_simple_parallel/fast:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/build
.PHONY : sdd_simple_parallel/fast

#=============================================================================
# Target rules for targets named strip_comments

# Build rule for target.
strip_comments: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 strip_comments
.PHONY : strip_comments

# fast build rule for target.
strip_comments/fast:
	$(MAKE) -f CMakeFiles/strip_comments.dir/build.make CMakeFiles/strip_comments.dir/build
.PHONY : strip_comments/fast

build_source_target_maps.o: build_source_target_maps.cc.o
.PHONY : build_source_target_maps.o

# target to build an object file
build_source_target_maps.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/build_source_target_maps.cc.o
.PHONY : build_source_target_maps.cc.o

build_source_target_maps.i: build_source_target_maps.cc.i
.PHONY : build_source_target_maps.i

# target to preprocess a source file
build_source_target_maps.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/build_source_target_maps.cc.i
.PHONY : build_source_target_maps.cc.i

build_source_target_maps.s: build_source_target_maps.cc.s
.PHONY : build_source_target_maps.s

# target to generate assembly for a file
build_source_target_maps.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/build_source_target_maps.cc.s
.PHONY : build_source_target_maps.cc.s

convergence_rates.o: convergence_rates.cc.o
.PHONY : convergence_rates.o

# target to build an object file
convergence_rates.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/convergence_rates.cc.o
.PHONY : convergence_rates.cc.o

convergence_rates.i: convergence_rates.cc.i
.PHONY : convergence_rates.i

# target to preprocess a source file
convergence_rates.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/convergence_rates.cc.i
.PHONY : convergence_rates.cc.i

convergence_rates.s: convergence_rates.cc.s
.PHONY : convergence_rates.s

# target to generate assembly for a file
convergence_rates.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/convergence_rates.cc.s
.PHONY : convergence_rates.cc.s

darcy_K_inverse.o: darcy_K_inverse.cc.o
.PHONY : darcy_K_inverse.o

# target to build an object file
darcy_K_inverse.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/darcy_K_inverse.cc.o
.PHONY : darcy_K_inverse.cc.o

darcy_K_inverse.i: darcy_K_inverse.cc.i
.PHONY : darcy_K_inverse.i

# target to preprocess a source file
darcy_K_inverse.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/darcy_K_inverse.cc.i
.PHONY : darcy_K_inverse.cc.i

darcy_K_inverse.s: darcy_K_inverse.cc.s
.PHONY : darcy_K_inverse.s

# target to generate assembly for a file
darcy_K_inverse.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/darcy_K_inverse.cc.s
.PHONY : darcy_K_inverse.cc.s

darcy_exact_solution_2.o: darcy_exact_solution_2.cc.o
.PHONY : darcy_exact_solution_2.o

# target to build an object file
darcy_exact_solution_2.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/darcy_exact_solution_2.cc.o
.PHONY : darcy_exact_solution_2.cc.o

darcy_exact_solution_2.i: darcy_exact_solution_2.cc.i
.PHONY : darcy_exact_solution_2.i

# target to preprocess a source file
darcy_exact_solution_2.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/darcy_exact_solution_2.cc.i
.PHONY : darcy_exact_solution_2.cc.i

darcy_exact_solution_2.s: darcy_exact_solution_2.cc.s
.PHONY : darcy_exact_solution_2.s

# target to generate assembly for a file
darcy_exact_solution_2.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/darcy_exact_solution_2.cc.s
.PHONY : darcy_exact_solution_2.cc.s

flux_function_2.o: flux_function_2.cc.o
.PHONY : flux_function_2.o

# target to build an object file
flux_function_2.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/flux_function_2.cc.o
.PHONY : flux_function_2.cc.o

flux_function_2.i: flux_function_2.cc.i
.PHONY : flux_function_2.i

# target to preprocess a source file
flux_function_2.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/flux_function_2.cc.i
.PHONY : flux_function_2.cc.i

flux_function_2.s: flux_function_2.cc.s
.PHONY : flux_function_2.s

# target to generate assembly for a file
flux_function_2.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/flux_function_2.cc.s
.PHONY : flux_function_2.cc.s

interfacial_pressure_function_2.o: interfacial_pressure_function_2.cc.o
.PHONY : interfacial_pressure_function_2.o

# target to build an object file
interfacial_pressure_function_2.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/interfacial_pressure_function_2.cc.o
.PHONY : interfacial_pressure_function_2.cc.o

interfacial_pressure_function_2.i: interfacial_pressure_function_2.cc.i
.PHONY : interfacial_pressure_function_2.i

# target to preprocess a source file
interfacial_pressure_function_2.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/interfacial_pressure_function_2.cc.i
.PHONY : interfacial_pressure_function_2.cc.i

interfacial_pressure_function_2.s: interfacial_pressure_function_2.cc.s
.PHONY : interfacial_pressure_function_2.s

# target to generate assembly for a file
interfacial_pressure_function_2.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/interfacial_pressure_function_2.cc.s
.PHONY : interfacial_pressure_function_2.cc.s

map_linker.o: map_linker.cc.o
.PHONY : map_linker.o

# target to build an object file
map_linker.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/map_linker.cc.o
.PHONY : map_linker.cc.o

map_linker.i: map_linker.cc.i
.PHONY : map_linker.i

# target to preprocess a source file
map_linker.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/map_linker.cc.i
.PHONY : map_linker.cc.i

map_linker.s: map_linker.cc.s
.PHONY : map_linker.s

# target to generate assembly for a file
map_linker.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/map_linker.cc.s
.PHONY : map_linker.cc.s

mesh_data.o: mesh_data.cc.o
.PHONY : mesh_data.o

# target to build an object file
mesh_data.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/mesh_data.cc.o
.PHONY : mesh_data.cc.o

mesh_data.i: mesh_data.cc.i
.PHONY : mesh_data.i

# target to preprocess a source file
mesh_data.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/mesh_data.cc.i
.PHONY : mesh_data.cc.i

mesh_data.s: mesh_data.cc.s
.PHONY : mesh_data.s

# target to generate assembly for a file
mesh_data.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/mesh_data.cc.s
.PHONY : mesh_data.cc.s

my_l2_norm.o: my_l2_norm.cc.o
.PHONY : my_l2_norm.o

# target to build an object file
my_l2_norm.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/my_l2_norm.cc.o
.PHONY : my_l2_norm.cc.o

my_l2_norm.i: my_l2_norm.cc.i
.PHONY : my_l2_norm.i

# target to preprocess a source file
my_l2_norm.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/my_l2_norm.cc.i
.PHONY : my_l2_norm.cc.i

my_l2_norm.s: my_l2_norm.cc.s
.PHONY : my_l2_norm.s

# target to generate assembly for a file
my_l2_norm.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/my_l2_norm.cc.s
.PHONY : my_l2_norm.cc.s

parallel_map_linker.o: parallel_map_linker.cc.o
.PHONY : parallel_map_linker.o

# target to build an object file
parallel_map_linker.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/parallel_map_linker.cc.o
.PHONY : parallel_map_linker.cc.o

parallel_map_linker.i: parallel_map_linker.cc.i
.PHONY : parallel_map_linker.i

# target to preprocess a source file
parallel_map_linker.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/parallel_map_linker.cc.i
.PHONY : parallel_map_linker.cc.i

parallel_map_linker.s: parallel_map_linker.cc.s
.PHONY : parallel_map_linker.s

# target to generate assembly for a file
parallel_map_linker.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/parallel_map_linker.cc.s
.PHONY : parallel_map_linker.cc.s

plot_i_data.o: plot_i_data.cc.o
.PHONY : plot_i_data.o

# target to build an object file
plot_i_data.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/plot_i_data.cc.o
.PHONY : plot_i_data.cc.o

plot_i_data.i: plot_i_data.cc.i
.PHONY : plot_i_data.i

# target to preprocess a source file
plot_i_data.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/plot_i_data.cc.i
.PHONY : plot_i_data.cc.i

plot_i_data.s: plot_i_data.cc.s
.PHONY : plot_i_data.s

# target to generate assembly for a file
plot_i_data.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/plot_i_data.cc.s
.PHONY : plot_i_data.cc.s

sdd_simple_parallel.o: sdd_simple_parallel.cc.o
.PHONY : sdd_simple_parallel.o

# target to build an object file
sdd_simple_parallel.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/sdd_simple_parallel.cc.o
.PHONY : sdd_simple_parallel.cc.o

sdd_simple_parallel.i: sdd_simple_parallel.cc.i
.PHONY : sdd_simple_parallel.i

# target to preprocess a source file
sdd_simple_parallel.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/sdd_simple_parallel.cc.i
.PHONY : sdd_simple_parallel.cc.i

sdd_simple_parallel.s: sdd_simple_parallel.cc.s
.PHONY : sdd_simple_parallel.s

# target to generate assembly for a file
sdd_simple_parallel.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/sdd_simple_parallel.cc.s
.PHONY : sdd_simple_parallel.cc.s

stokes_exact_solution_2.o: stokes_exact_solution_2.cc.o
.PHONY : stokes_exact_solution_2.o

# target to build an object file
stokes_exact_solution_2.cc.o:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/stokes_exact_solution_2.cc.o
.PHONY : stokes_exact_solution_2.cc.o

stokes_exact_solution_2.i: stokes_exact_solution_2.cc.i
.PHONY : stokes_exact_solution_2.i

# target to preprocess a source file
stokes_exact_solution_2.cc.i:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/stokes_exact_solution_2.cc.i
.PHONY : stokes_exact_solution_2.cc.i

stokes_exact_solution_2.s: stokes_exact_solution_2.cc.s
.PHONY : stokes_exact_solution_2.s

# target to generate assembly for a file
stokes_exact_solution_2.cc.s:
	$(MAKE) -f CMakeFiles/sdd_simple_parallel.dir/build.make CMakeFiles/sdd_simple_parallel.dir/stokes_exact_solution_2.cc.s
.PHONY : stokes_exact_solution_2.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... debug"
	@echo "... distclean"
	@echo "... edit_cache"
	@echo "... info"
	@echo "... rebuild_cache"
	@echo "... release"
	@echo "... run"
	@echo "... runclean"
	@echo "... sdd_simple_parallel"
	@echo "... strip_comments"
	@echo "... build_source_target_maps.o"
	@echo "... build_source_target_maps.i"
	@echo "... build_source_target_maps.s"
	@echo "... convergence_rates.o"
	@echo "... convergence_rates.i"
	@echo "... convergence_rates.s"
	@echo "... darcy_K_inverse.o"
	@echo "... darcy_K_inverse.i"
	@echo "... darcy_K_inverse.s"
	@echo "... darcy_exact_solution_2.o"
	@echo "... darcy_exact_solution_2.i"
	@echo "... darcy_exact_solution_2.s"
	@echo "... flux_function_2.o"
	@echo "... flux_function_2.i"
	@echo "... flux_function_2.s"
	@echo "... interfacial_pressure_function_2.o"
	@echo "... interfacial_pressure_function_2.i"
	@echo "... interfacial_pressure_function_2.s"
	@echo "... map_linker.o"
	@echo "... map_linker.i"
	@echo "... map_linker.s"
	@echo "... mesh_data.o"
	@echo "... mesh_data.i"
	@echo "... mesh_data.s"
	@echo "... my_l2_norm.o"
	@echo "... my_l2_norm.i"
	@echo "... my_l2_norm.s"
	@echo "... parallel_map_linker.o"
	@echo "... parallel_map_linker.i"
	@echo "... parallel_map_linker.s"
	@echo "... plot_i_data.o"
	@echo "... plot_i_data.i"
	@echo "... plot_i_data.s"
	@echo "... sdd_simple_parallel.o"
	@echo "... sdd_simple_parallel.i"
	@echo "... sdd_simple_parallel.s"
	@echo "... stokes_exact_solution_2.o"
	@echo "... stokes_exact_solution_2.i"
	@echo "... stokes_exact_solution_2.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

