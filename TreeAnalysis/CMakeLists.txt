cmake_minimum_required(VERSION 2.6)

project(MVALeptonID) 

#options here
set(CMAKE_CXX_FLAGS "-Wall -march=native")

#get ROOT stuff here
exec_program(root-config ARGS "--incdir" OUTPUT_VARIABLE ROOT_INCLUDE_DIR)
include_directories(${ROOT_INCLUDE_DIR})

exec_program(root-config ARGS "--libs" OUTPUT_VARIABLE ROOT_LIBS)
set(ROOT_LIBS ${ROOT_LIBS} "-lRooFit -lRooFitCore -lMinuit -lFoam")
message(STATUS ${ROOT_LIBS})

set(DataFormats ../DataFormats)
set(DataFormats_interface ${DataFormats}/interface)
set(DataFormats_src ${DataFormats}/src)
set(TreeAnalysis_interface interface)
set(BIN bin)
set(SRC src)

set(CMAKE_INCLUDE_PATH ${DataFormats_interface} ${TreeAnalysis_interface}) 
set(EventDF_HEADERS ${DataFormats_interface}/*.h)
set(EventDF_LINKDEF ${DataFormats_src}/STDToolsLinkDef.h)
set(EventDF_DICTIONARY STDTools.cxx) 

include(FindROOT.cmake)
ROOT_GENERATE_DICTIONARY("${EventDF_HEADERS}" "${EventDF_LINKDEF}" "${EventDF_DICTIONARY}" "${include_directories}")

include_directories(${CMAKE_INCLUDE_PATH})
set(EventAnalyzer_SRCS ${EventDF_DICTIONARY} ${SRC}/EventAnalyzer.cc ${SRC}/ZZWAnalyzer.cc)
add_executable(${BIN}/eventAnalyzer ${SRC}/eventAnalyzer.cpp ${EventAnalyzer_SRCS})
target_link_libraries(bin/eventAnalyzer ${ROOT_LIBS})






