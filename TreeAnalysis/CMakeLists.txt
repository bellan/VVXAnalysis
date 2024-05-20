##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################

cmake_minimum_required(VERSION 3.4)

project(VVXAnalysis) 

set(CMAKE_MACOSX_RPATH 0)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---locate Python lib (used to embed python scripts in C)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

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
set(EventDF_LINKDEF ${DataFormats_src}/STDToolsLinkDef.h)

include_directories(${CMAKE_INCLUDE_PATH} ${DataFormats_interface} ${BASE})

set(EventDF_DICTIONARY G__TreeAnalysis)
ROOT_GENERATE_DICTIONARY(${EventDF_DICTIONARY} ${EventDF_HEADERS} LINKDEF "${EventDF_LINKDEF}")


set(EventAnalyzer_SRCS ${EventDF_DICTIONARY} ${SRC}/SampleInfo.cc ${SRC}/EventAnalyzer.cc ${Commons_src}/SignalDefinitions.cc ${DataFormats_src}/RegionTypes.cc ${Commons_src}/StringTools.cc ${Commons_src}/GenVBHelper.cc)


set(Analyzers_SRCS)
file(GLOB Analyzers_SRCS ${SRC}/*.cc)
set(VVXAnalyzer_SRCS ${EventAnalyzer_SRCS} ${Analyzers_SRCS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

add_library(TreeAnalysis SHARED ${EventDF_SOURCES} G__TreeAnalysis.cxx)
target_link_libraries(TreeAnalysis ${ROOT_LIBRARIES})
add_executable(eventAnalyzer ${BIN}/eventAnalyzer.cpp ${SRC}/AnalysisFactory.cc ${VVXAnalyzer_SRCS})
target_link_libraries(eventAnalyzer TreeAnalysis)
target_link_libraries(eventAnalyzer ${PYTHON_LIBRARIES})

add_executable(treeFilter ${BIN}/treeFilter.cpp ${EventDF_DICTIONARY})
target_link_libraries(treeFilter TreeAnalysis)

add_executable(genEventAnalyzer ${BIN}/genEventAnalyzer.cpp ${SRC}/GenEventAnalyzer.cc ${EventAnalyzer_SRCS})
target_link_libraries(genEventAnalyzer TreeAnalysis)
