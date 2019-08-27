# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cristiano/VVXAnalysis/TreeAnalysis

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cristiano/VVXAnalysis/TreeAnalysis

# Utility rule file for G__TreeAnalysis.

# Include the progress variables for this target.
include CMakeFiles/G__TreeAnalysis.dir/progress.make

CMakeFiles/G__TreeAnalysis: G__TreeAnalysis.cxx
CMakeFiles/G__TreeAnalysis: libTreeAnalysis_rdict.pcm
CMakeFiles/G__TreeAnalysis: libTreeAnalysis.rootmap


G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/src/STDToolsLinkDef.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/DiBoson.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/TypeDefs.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/GenEventWeights.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/MELA.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Jet.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Boson.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/GenStatusBit.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Lepton.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Electron.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Particle.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/DiBoson.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/TypeDefs.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/GenEventWeights.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/MELA.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Jet.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Boson.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/GenStatusBit.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Lepton.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Electron.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/interface/Particle.h
G__TreeAnalysis.cxx: /home/cristiano/VVXAnalysis/DataFormats/src/STDToolsLinkDef.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/cristiano/VVXAnalysis/TreeAnalysis/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating G__TreeAnalysis.cxx, libTreeAnalysis_rdict.pcm, libTreeAnalysis.rootmap"
	/usr/bin/cmake -E env LD_LIBRARY_PATH=/usr/local/root/lib:/usr/local/root/lib /usr/local/root/bin/rootcling -v2 -f G__TreeAnalysis.cxx -s /home/cristiano/VVXAnalysis/TreeAnalysis/libTreeAnalysis.so -rml libTreeAnalysis.so -rmf /home/cristiano/VVXAnalysis/TreeAnalysis/libTreeAnalysis.rootmap -I/usr/local/root/include -I/home/cristiano/VVXAnalysis/TreeAnalysis -I/home/cristiano/VVXAnalysis/TreeAnalysis/../.. -I/home/cristiano/VVXAnalysis/TreeAnalysis/interface -I/opt/local/include -I/home/cristiano/VVXAnalysis/TreeAnalysis/../DataFormats/interface /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/DiBoson.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/TypeDefs.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/GenEventWeights.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/MELA.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Jet.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Boson.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/GenStatusBit.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Lepton.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Electron.h /home/cristiano/VVXAnalysis/TreeAnalysis/../..//VVXAnalysis/DataFormats/interface/Particle.h /home/cristiano/VVXAnalysis/TreeAnalysis/../DataFormats/src/STDToolsLinkDef.h

libTreeAnalysis_rdict.pcm: G__TreeAnalysis.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libTreeAnalysis_rdict.pcm

libTreeAnalysis.rootmap: G__TreeAnalysis.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libTreeAnalysis.rootmap

G__TreeAnalysis: CMakeFiles/G__TreeAnalysis
G__TreeAnalysis: G__TreeAnalysis.cxx
G__TreeAnalysis: libTreeAnalysis_rdict.pcm
G__TreeAnalysis: libTreeAnalysis.rootmap
G__TreeAnalysis: CMakeFiles/G__TreeAnalysis.dir/build.make

.PHONY : G__TreeAnalysis

# Rule to build all files generated by this target.
CMakeFiles/G__TreeAnalysis.dir/build: G__TreeAnalysis

.PHONY : CMakeFiles/G__TreeAnalysis.dir/build

CMakeFiles/G__TreeAnalysis.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/G__TreeAnalysis.dir/cmake_clean.cmake
.PHONY : CMakeFiles/G__TreeAnalysis.dir/clean

CMakeFiles/G__TreeAnalysis.dir/depend:
	cd /home/cristiano/VVXAnalysis/TreeAnalysis && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cristiano/VVXAnalysis/TreeAnalysis /home/cristiano/VVXAnalysis/TreeAnalysis /home/cristiano/VVXAnalysis/TreeAnalysis /home/cristiano/VVXAnalysis/TreeAnalysis /home/cristiano/VVXAnalysis/TreeAnalysis/CMakeFiles/G__TreeAnalysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/G__TreeAnalysis.dir/depend

