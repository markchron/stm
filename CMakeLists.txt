cmake_minimum_required(VERSION 2.8.4)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug
		CACHE STRING
		"Choose the build type from Debug Release profiles"
		FORCE
		)
endif()
set(CMAKE_INSTALL_PREFIX "~/local"
	CACHE STRING
	"The install prefix (default ~/local)."
	)
option(CMAKE_VERBOSE_MAKEFILE "make verbose off | make VERBOSE=1" OFF)
## switch on|off Intel
## cmake -DST_WITH_INTEL:BOOL=ON path/to/CMakeLists.txt
option(ST_WITH_INTEL "compile by Intel compiler." ON)
option(ST_WITH_OPENMP "active OpenMP parallelization." ON)
include(${CMAKE_SOURCE_DIR}/cmake/options.cmake)

project(stm1023 Fortran)
message(STATUS "Build type 		: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install path 	: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMake path		: ${CMAKE_MODULE_PATH}")
#message(STATUS "Project source 	: ${CMAKE_SOURCE_DIR} ") # = ${PROJECT_SOURCE_DIR}
#message(STATUS "Project build	: ${PROJECT_BINARY_DIR}")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
include(conflags)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
