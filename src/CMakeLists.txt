## add all the source file in current folder to the SRC_LIST
## empty file lists
aux_source_directory(./ SRC_LIST)
message(STATUS "source files SRC_LIST : ${SRC_LIST}")

## add all the files started with st to the ST_SRC_LIST
# file(GLOB ST_SRC_LIST "st*") 
## add support the recursive searching file in subfolders
file(GLOB_RECURSE ST_SRC_LIST "./st*")
## we DONOT recommend using GLOB to collect a list of source files from your
## source tree. If no CMakeList.txt file changes when a source is added or
## removed then the generated build system cannot know when to ask CMake to
## regenerate. 
## MUST RUN cmake again after adding or removing files

message(STATUS "source files ST_SRC_LIST: ${ST_SRC_LIST}")
message(STATUS "libraries  ST_LIBS_LIST : ${ST_LIBS_LIST}")

## on wellbore
set(WLL_SRCS 
	${PROJECT_SOURCE_DIR}/src/sthead.f90
	${PROJECT_SOURCE_DIR}/src/stgeom.f90
	${PROJECT_SOURCE_DIR}/src/stwell.f90
	${PROJECT_SOURCE_DIR}/src/testwell.f90
	)
add_executable(TestWll ${WLL_SRCS})
target_link_libraries(TestWll ${ST_LIBS_LIST})

add_executable(stm ${ST_SRC_LIST})
target_link_libraries(stm ${ST_LIBS_LIST})
## create a folder "executables" and adds target
## project (stm.vcproj) under it
set_property(TARGET stm PROPERTY FOLDER "executables")

## adds logic to INSTALL.vcproj to copy stm.exe to destination direcory
install(TARGETS stm
	RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin)
