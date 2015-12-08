	  module stmcomapi
      use stmdatpol, ONLY : MASTER
	  integer :: err
      include 'mpif.h'
      
      contains
      subroutine pmpi_fname(infile, lth, trunc)
      integer, intent(inout)			:: lth
      character(len=lth), intent(inout)	:: infile
      integer, intent(inout) 			:: trunc

      call MPI_BCAST(lth, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, err)
      call MPI_BCAST(infile, lth, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, err)
      
      call MPI_BCAST(trunc, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, err)
	  end subroutine pmpi_fname

	  end module stmcomapi

