cmake_minimum_required(VERSION 2.8.4)
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

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
option(ST_WITH_MPI "active OPENMPI parallelization." ON)
option(ST_WITH_OPENMP "active OpenMP parallelization." ON)
include(${CMAKE_SOURCE_DIR}/cmake/options.cmake)

message(STATUS "CMake module	: ${CMAKE_MODULE_PATH}")
message(STATUS "CMake prefix	: ${CMAKE_PREFIX_PATH}")
#message(STATUS "PATH 	: $ENV{PATH}")
#enable_language(Fortran)

project(stm1023 Fortran)
message(STATUS "Build type 		: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install path 	: ${CMAKE_INSTALL_PREFIX}")
#message(STATUS "Project source 	: ${CMAKE_SOURCE_DIR} ") # = ${PROJECT_SOURCE_DIR}
#message(STATUS "Project build	: ${PROJECT_BINARY_DIR}")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

include(conflags)
include(confprefix)
if(ST_WITH_MPI) 	
	find_package(MPI REQUIRED)
	include_directories(${MPI_Fortran_INCLUDE_PATH})
	link_directories(${MPI_Fortran_LIBRARIES})
	set(ST_LIBS_LIST ${MPI_Fortran_LIBRARIES})
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(test)

#find_package(LATEX PATHS ${ST_LATEX_DIR})  ## latex.conf
find_package(LATEX)
if(UNIX AND LATEX_FOUND)
	add_subdirectory(doc)
elseif(NOT LATEX_FOUND)
	message(STATUS "MIKTEX_BINARY_PATH 	: ${MIKTEX_BINARY_PATH}")
endif() #latex
