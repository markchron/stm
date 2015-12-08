file(GLOB test_dat_files
	./*.dat 
	)
file(COPY ${test_dat_files}
	DESTINATION ${PROJECT_BINARY_DIR}/test
	)

macro(do_mpi_test arg np)
	add_test(NAME MP-${np}-${arg}
		COMMAND ${MPIEXEC} -n ${np} ${PROJECT_BINARY_DIR}/bin/stm ${arg}
		)
endmacro (do_mpi_test)

add_test(NAME spe4b
	COMMAND ${PROJECT_BINARY_DIR}/bin/stm spe4b.dat
	)

do_mpi_test(spe4b.dat 2)
