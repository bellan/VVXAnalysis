##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################

cmake_minimum_required(VERSION 2.6)

project(VVXAnalysis) 

#options here
set(CMAKE_CXX_FLAGS "-Wall -march=native -std=c++0x")

#get ROOT stuff here
exec_program(root-config ARGS "--incdir" OUTPUT_VARIABLE ROOT_INCLUDE_DIR)
include_directories(${ROOT_INCLUDE_DIR})

exec_program(root-config ARGS "--libs" OUTPUT_VARIABLE ROOT_LIBS)
set(ROOT_LIBS ${ROOT_LIBS} "-lRooFit -lRooFitCore -lMinuit -lFoam")
message(STATUS ${ROOT_LIBS})

set(DataFormats ../DataFormats)
set(DataFormats_src ${DataFormats}/src)
set(Commons_src ../Commons/src)


set(BASE ../../)

set(TreeAnalysis_interface interface)
set(BIN bin)
set(SRC src)

set(CMAKE_INCLUDE_PATH ${BASE} ${TreeAnalysis_interface}) 
set(EventDF_HEADERS ${BASE}/VVXAnalysis/DataFormats/interface/*.h)
set(EventDF_LINKDEF ${DataFormats_src}/STDToolsLinkDef.h)
set(EventDF_DICTIONARY STDTools.cxx)

include(FindROOT.cmake)
ROOT_GENERATE_DICTIONARY("${EventDF_HEADERS}" "${EventDF_LINKDEF}" "${EventDF_DICTIONARY}" "${include_directories}")

include_directories(${CMAKE_INCLUDE_PATH})
set(EventAnalyzer_SRCS ${EventDF_DICTIONARY} ${SRC}/MCInfo.cc ${SRC}/EventAnalyzer.cc ${Commons_src}/SignalDefinitions.cc ${Commons_src}/LeptonScaleFactors.cc ${Commons_src}/RegionTypes.cc ${Commons_src}/StringTools.cc)

#Set here the location of the source code of your analysis
set(VVXAnalyzer_SRCS ${EventAnalyzer_SRCS} ${SRC}/VVXAnalyzer.cc ${SRC}/ZZWAnalyzer.cc ${SRC}/ZZSAnalyzer.cc ${SRC}/ZZjAnalyzer.cc ${SRC}/ZZWSRDefinition.cc ${SRC}/ZZMCAnalyzer.cc ${SRC}/ZZRecoAnalyzer.cc ${SRC}/FakeRateAnalyzer.cc ${SRC}/ZZjGenAnalyzer.cc ${SRC}/VBSMCAnalyzer.cc ${SRC}/VBSRecoAnalyzer.cc)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

add_executable(eventAnalyzer ${BIN}/eventAnalyzer.cpp ${SRC}/AnalysisFactory.cc ${VVXAnalyzer_SRCS})
target_link_libraries(eventAnalyzer ${ROOT_LIBS})

add_executable(treeFilter ${BIN}/treeFilter.cpp ${EventDF_DICTIONARY})
target_link_libraries(treeFilter ${ROOT_LIBS})

add_executable(genEventAnalyzer ${BIN}/genEventAnalyzer.cpp ${SRC}/GenEventAnalyzer.cc ${EventAnalyzer_SRCS})
target_link_libraries(genEventAnalyzer ${ROOT_LIBS})





