# Try to find the METIS libraries
# METIS_FOUND - system has METIS lib
# METIS_INCLUDE_DIR - the METIS include directory
# METIS_LIBRARIES - the METIS libraries

if(WIN32)
        ## todo
        #set(GLEW_LIBRARIES "glew32")
        #set(GLEW_FOUND CACHE INTERNAL TRUE)
endif(WIN32)

if(UNIX)
	find_path(METIS_INCLUDE_DIR 
		NAMES metis.h
		PATHS ${ST_METIS_DIR}
		/usr
		/opt
		/opt/local
		/sw
		${HOME}/usr
		PATH_SUFFIXES include
		)

	find_library(METIS_LIBRARIES 
		NAMES libmetis.a 							#metis libmetis, INTERFACE WT Fortran, must be static library
		PATHS ${ST_METIS_DIR}
		/usr
		/usr/local
		/opt
		/opt/local
		/sw
		${HOME}/usr
		PATH_SUFFIXES lib
		)
	
	if(METIS_INCLUDE_DIR AND METIS_LIBRARIES)
		set(METIS_FOUND TRUE)
		add_definitions(-DENABLE_METIS)
	endif()

	if(METIS_FOUND)
		if(NOT METIS_FIND_QUIETLY)
			message(STATUS "Found METIS : ${METIS_LIBRARIES} ${METIS_INCLUDE_DIR}")
		endif()
	else()
		if(METIS_FIND_REQUIRED)
			message(FATAL_ERROR "Could not find METIS library. ${METIS_LIBRARIES} ${METIS_INCLUDE_DIR}")
        endif()
	endif()
endif(UNIX)
