# CMakeList.txt: splitting_filter 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# 将源代码添加到此项目的可执行文件。
aux_source_directory(. DIR_SRCS)
add_executable(splitting_filter ${DIR_SRCS})	

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/DEBUG/)	# 设置可执行文件输出目录
##SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})	# 设置可执行文件输出目录
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${PROJECT_SOURCE_DIR})
#include_directories (${PROJECT_SOUREC_DIR}/include)
#install(FILES ${PROJECT_SOURCE_DIR}*.wav DESTINATION ${PROJECT_SOURCE_DIR}DEBUG/)