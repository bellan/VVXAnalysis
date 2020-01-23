##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################

cmake_minimum_required(VERSION 2.6)

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
set(EventDF_LINKDEF ${DataFormats_src}/STDToolsLinkDef.h)

include_directories(${CMAKE_INCLUDE_PATH} ${DataFormats_interface} ${BASE})

ROOT_GENERATE_DICTIONARY(G__TreeAnalysis "${EventDF_HEADERS}" LINKDEF "${EventDF_LINKDEF}")



set(EventAnalyzer_SRCS ${EventDF_DICTIONARY} ${SRC}/MCInfo.cc ${SRC}/EventAnalyzer.cc ${Commons_src}/SignalDefinitions.cc ${Commons_src}/LeptonScaleFactors.cc ${Commons_src}/RegionTypes.cc ${Commons_src}/StringTools.cc)

#Set here the location of the source code of your analysis
set(VVXAnalyzer_SRCS ${EventAnalyzer_SRCS} ${SRC}/VVXAnalyzer.cc ${SRC}/ZZWAnalyzer.cc ${SRC}/WZAnalyzer.cc ${SRC}/ZZSAnalyzer.cc ${SRC}/ZZjAnalyzer.cc ${SRC}/ZZWSRDefinition.cc ${SRC}/ZZMCAnalyzer.cc ${SRC}/ZZRecoAnalyzer.cc ${SRC}/FakeRateAnalyzer.cc ${SRC}/ZZjGenAnalyzer.cc ${SRC}/VBSMCAnalyzer.cc ${SRC}/VBSRecoAnalyzer.cc ${SRC}/WlllnuAnalyzer.cc ${SRC}/ZVAnalyzer.cc ${SRC}/WWosAnalyzer.cc ${SRC}/WZZAnalyzer.cc ${SRC}/DBGenerator.cc)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

add_library(TreeAnalysis SHARED G__TreeAnalysis.cxx)
target_link_libraries(TreeAnalysis ${ROOT_LIBRARIES})
add_executable(eventAnalyzer ${BIN}/eventAnalyzer.cpp ${SRC}/AnalysisFactory.cc ${VVXAnalyzer_SRCS})
target_link_libraries(eventAnalyzer TreeAnalysis)

add_executable(treeFilter ${BIN}/treeFilter.cpp ${EventDF_DICTIONARY})
target_link_libraries(treeFilter TreeAnalysis)

add_executable(genEventAnalyzer ${BIN}/genEventAnalyzer.cpp ${SRC}/GenEventAnalyzer.cc ${EventAnalyzer_SRCS})
target_link_libraries(genEventAnalyzer TreeAnalysis)





