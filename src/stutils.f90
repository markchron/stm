	  module stmutils
	  use stmheader
	  use stmdatpol
	  use stmdatnpt
	  use stmcomapi
	  use stmoutput
	  
      contains
	  subroutine st_init_memo
	  include 'mpif.h'
      integer :: mprov
	  nprocs = 1
	  rank = 0
	  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mprov, stErrs(6))
	  if(stErrs(6) == 0) then
		  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, stErrs(7))
		  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, stErrs(8))
	  endif
      deck = rank + 1

      call omp_system_optimal_num_threads
	  ! input file name
	  call get_ifile
	  if(nprocs > 1 ) call pmpi_fname(stifnm, stiflen, stiftrunc)
	  ! create the *.out files & *.log
	  call dprt_create

	  ! check input file exists or not 
	  call serial_verify_ifile

	  ! initialize memory to store controling variables
      call datpol_init_scalar
	  ! define the controlling variables 
	  call npt_init
      ! asks for memory
      call datpol_init_memo
	  end subroutine st_init_memo

	  subroutine st_reservoir_init
      ! read in field & wellbore info.
      call npt_set_properties
      ! update vector or control variables
	  call datpol_update_index

	  end subroutine st_reservoir_init
! PURPOSE:
! distribute the problem between processors by distribute grid cells
      subroutine st_dist
      call set_partition_dist
      ! master processor collects the local information
      call pmpi_master_collect_i(stIcsv(110:115), 6, nprocs,  pmInmas)
      end subroutine st_dist
! PURPOSE:
! terminate the parallel program
! close the opening files
! release the dynamic memory
	  subroutine st_release_memo
      use stmoutput, only : dprt_errmsg
      include 'mpif.h'
	  close(I_UNIT_5, IOSTAT = stErrs(9))
	  close(FUNIT_OUT, IOSTAT = stErrs(10))
       ! print the error msgs
	  call dprt_errmsg
	  call datpol_free
	  call MPI_FINALIZE(stErrs(6)) 
      end subroutine st_release_memo
!----------------------------------------------------------------------      
! OPENMP
!----------------------------------------------------------------------
! PURPOSE:
! set up multi-threads number
      subroutine omp_system_optimal_num_threads
      use omp_lib
      nthread = omp_get_max_threads()
      call OMP_SET_NUM_THREADS(nthread)
      end subroutine omp_system_optimal_num_threads
	  end module stmutils
