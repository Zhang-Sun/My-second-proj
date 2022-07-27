PROJECT(AdjLib)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

file(GLOB SRC_FILES "*.cpp")
add_library(AdjLib ${SRC_FILES})
target_link_libraries(AdjLib LogLib CmnLib)
