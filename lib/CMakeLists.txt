if(NOT METIS_FOUND)
	message(STATUS "Compile static METIS library with supplied source files.")
	add_subdirectory(metis-5.1.0)
	#list(APPEND METIS_LIB "metis")
	#	message(STATUS "METIS_LIB ${METIS_LIB}")
endif()
