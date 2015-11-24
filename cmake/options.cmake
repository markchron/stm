if(ST_WITH_INTEL)
	find_program(CMAKE_C_COMPILER NAMES icc)
	find_program(CMAKE_CXX_COMPILER NAMES icpc)
	find_program(CMAKE_Fortran_COMPILER NAMES ifort)
	find_program(CMAKE_AR NAMES xiar)
	find_program(CMAKE_LINKER NAMES xild)
	if(CMAKE_C_COMPILER MATCHES CMAKE_C_COMPILER-NOTFOUND OR
			CMAKE_CXX_COMPILER MATCHES CMAKE_CXX_COMPILER-NOTFOUND OR
			CMAKE_Fortran_COMPILER MATCHES CMAKE_Fortran_COMPILER-NOTFOUND OR
			CMAKE_AR MATCHES CMAKE_AR-NOTFOUND OR
			CMAKE_LINKER MATCHES CMAKE_LINKER-NOTFOUND)
		message(WARNING "Intel ${ST_WITH_INTEL}, failed to find its compiler at default PATH." )
		#		message(STATUS	"Leave CMake to find avaliable compilers.")
	endif()
endif()

