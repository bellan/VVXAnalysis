##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################

cmake_minimum_required(VERSION 3.4)

project(VVXAnalysis) 

set(CMAKE_MACOSX_RPATH 0)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


set(DataFormats ../DataFormats)
set(DataFormats_src ${DataFormats}/src)
set(DataFormats_interface ${DataFormats}/interface)
set(Commons_src ../Commons/src)

set(BASE ../../)

set(TreeAnalysis_interface interface)
set(BIN bin)
set(SRC src)

set(BOOST /opt/local/include/)

set(CMAKE_INCLUDE_PATH ${BASE} ${TreeAnalysis_interface} ${BOOST}) 
set(EventDF_HEADERS ${BASE}/VVXAnalysis/DataFormats/interface/*.h)
set(EventDF_SOURCES ${DataFormats_src}/Particle.cc ${DataFormats_src}/Lepton.cc ${DataFormats_src}/Jet.cc ${DataFormats_src}/Photon.cc ${DataFormats_src}/Electron.cc ${DataFormats_src}/Boson.cc ${DataFormats_src}/DiBoson.cc ${DataFormats_src}/RegionTypes.cc ${DataFormats_src}/loader.C)

include_directories(${CMAKE_INCLUDE_PATH} ${DataFormats_interface} ${BASE})

set(EventDF_DICTIONARY VVXAnalysisDataFormats)
REFLEX_GENERATE_DICTIONARY(${EventDF_DICTIONARY} ${DataFormats_src}/classes.h SELECTION ${DataFormats_src}/classes_def.xml)
set_property(TARGET ${EventDF_DICTIONARY}-dictgen APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${EventDF_DICTIONARY}_rdict.pcm)


set(EventAnalyzer_SRCS ${SRC}/SampleInfo.cc ${SRC}/EventAnalyzer.cc ${Commons_src}/SignalDefinitions.cc ${DataFormats_src}/RegionTypes.cc ${Commons_src}/StringTools.cc ${Commons_src}/GenVBHelper.cc)


set(Analyzers_SRCS)
file(GLOB Analyzers_SRCS ${SRC}/*.cc)
list(REMOVE_ITEM Analyzers_SRCS ${SRC}/VZZAnalyzer.cc)
list(REMOVE_ITEM Analyzers_SRCS ${SRC}/DBGenerator.cc)
set(VVXAnalyzer_SRCS ${EventAnalyzer_SRCS} ${Analyzers_SRCS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

add_library(TreeAnalysis SHARED ${EventDF_SOURCES} ${EventDF_DICTIONARY}.cxx)
target_link_libraries(TreeAnalysis ${ROOT_LIBRARIES})
add_executable(eventAnalyzer ${BIN}/eventAnalyzer.cpp ${SRC}/AnalysisFactory.cc ${VVXAnalyzer_SRCS})
target_link_libraries(eventAnalyzer TreeAnalysis)
target_link_libraries(eventAnalyzer ${PYTHON_LIBRARIES})

add_executable(treeFilter ${BIN}/treeFilter.cpp ${EventDF_DICTIONARY})
target_link_libraries(treeFilter TreeAnalysis)

add_executable(genEventAnalyzer ${BIN}/genEventAnalyzer.cpp ${SRC}/GenEventAnalyzer.cc ${EventAnalyzer_SRCS})
target_link_libraries(genEventAnalyzer TreeAnalysis)
