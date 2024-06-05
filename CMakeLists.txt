cmake_minimum_required(VERSION 3.0)
project(JAEAExp)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Analysis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/ChannelData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/source/DetectorSetting.cpp
)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Analysis.h
    ${CMAKE_CURRENT_SOURCE_DIR}/source/ChannelData.h
    ${CMAKE_CURRENT_SOURCE_DIR}/source/DetectorSetting.h
)

find_package(ROOT REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
root_generate_dictionary(G__JAEAExp ${HEADERS} LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/source/LinkDef.h)
add_library(JAEAExp SHARED ${SOURCES} G__JAEAExp.cxx)
target_include_directories(JAEAExp PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(JAEAExp PUBLIC ${ROOT_LIBRARIES})
