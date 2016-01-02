      program stmexe
      use stmutils 
      use stmoutput
      implicit none
      ! set up mpi, openmp env.
      ! parse the input file name
      ! define problem size, nx,ny,nz
      ! update intermedia control size variables
      ! allocate the memory and pointers
      call st_init_memo
      ! define the field, wellbore info.
      ! update control var. & vectors based on known vectors
      call st_reservoir_init
      ! domain partitioning
      call st_dist
      ! output into the *.out
      call dprt_datpol
      ! release memo, nullify pointers
      ! print out errors
      ! close mpi env.
      call st_release_memo
      end program stmexe
