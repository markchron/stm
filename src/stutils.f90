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
	  ! input file name
	  call get_ifile
	  if(nprocs > 1 ) call pmpi_fname(stifnm, stiflen, stiftrunc)
	  ! create the *.out files 
	  call dprt_create
	  ! check input file exists or not 
	  call serial_verify_ifile

	  ! initialize memory to store controling variables
      call datpol_init_scalar
	  ! define the controlling variables 
	  call npt_init
	  end subroutine st_init_memo

	  subroutine st_reservoir_init
      call npt_set_properties

	  call datpol_update_index
      call dprt_csr(Ngcll, Ngedges, ivadj, iadjncy, detrans, "transmissibility", 1)

	  end subroutine st_reservoir_init
! PURPOSE:
! terminate the parallel program
! close the opening files
! release the dynamic memory
	  subroutine st_release_memo
      include 'mpif.h'
	  close(I_UNIT_5, IOSTAT = stErrs(9))
	  close(FUNIT_OUT, IOSTAT = stErrs(10))

	  if(rank == MASTER) call errmsg
	  call datpol_free
	  call MPI_FINALIZE(stErrs(6)) 
      end subroutine st_release_memo
	  end module stmutils
