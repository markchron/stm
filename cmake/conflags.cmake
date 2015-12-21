if(CMAKE_Fortran_COMPILER_WORKS)
	message(STATUS "Fortran compiler 	: ${CMAKE_SYSTEM} : ${CMAKE_Fortran_COMPILER_ID}")

	if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU") ## -Wtabs
		set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -ggdb -fbounds-check -fsignaling-nans -fno-f2c -ffpe-trap=zero,invalid,overflow")
		set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -pipe -funroll-all-loops")
		if(ST_WITH_OPENMP)
			set(ST_OPENMP_OPTS "-fopenmp")
		endif()
	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
		if(${UNIX})
			#			set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -check all -warn all -fpe-all=3")
			## attempt to use pointer when it is not associated with a target
			set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -check bounds -check uninit -check format -warn all -fpe-all=3")
			set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -vec-guard-write -fpconstant -extend-source -funroll-loops -align all -ip")
			if(ST_WITH_OPENMP)
				set(ST_OPENMP_OPTS "-openmp")
			endif()
		elseif(${WIN32})
			set(CMAKE_Fortran_FLAGS_DEBUG "/FA /debug:all /traceback /4Yb /warn:all /fpe-all:3")
			set(CMAKE_Fortran_FLAGS_RELEASE "/O3 /fast /QxHost /Qvec-guard-write -Qunroll /align /Qip")
			if(ST_WITH_OPENMP)
				set(ST_OPENMP_OPTS "/Qopenmp")
		    endif()		
		endif()
	elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "MSVC")
		message(SEND_ERROR "To do Fortran compiler flags for Visual Studio.")
	elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Clang")
		message(SEND_ERROR "To do Fortran compiler flags for APPLE CLANG.")
	else()
		if(${WIN32})
			set(CMAKE_Fortran_FLAGS_DEBUG "/FA /debug:all /traceback /4Yb /warn:all /fpe-all:3")
			set(CMAKE_Fortran_FLAGS_RELEASE "/O3 /fast /QxHost /Qvec-guard-write -Qunroll /align /Qip")
			if(ST_WITH_OPENMP)
				set(ST_OPENMP_OPTS "/Qopenmp")
		    endif()		
		else()
			message(SEND_ERROR "no-config Fortran compiler ${CMAKE_Fortran_COMPILER}")
		endif()
	endif()

	#	message(STATUS "Fortran flags 	: ${CMAKE_Fortran_FLAGS}")
	if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
		set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${ST_OPENMP_OPTS}")
		set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
		#list(APPEND CMAKE_Fortran_FLAGS_RELEASE ${ST_OPENMP_OPTS}) ## failed config
		set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${ST_OPENMP_OPTS}")
		set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
	else()
		message(SEND_ERROR "unknown build type, unknown compiler flags")
	endif()
	message(STATUS "compiler 		: ${CMAKE_Fortran_COMPILER}")
	message(STATUS "Fortran flags 	: ${CMAKE_Fortran_FLAGS}")
	
endif() ## CMAKE_Fortran_COMPILER_WORKS

