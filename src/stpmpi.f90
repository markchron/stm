      module stmcomapi
      use stmdatpol, ONLY : MASTER
      include 'mpif.h'
      integer :: err
      integer, dimension(:), allocatable        :: pmIsscnt, pmIssdis
      ! scatter operation: master processor distributes the property
      ! pmIsscnt, integer, Scatterv Sendout CouNTs
      !     * specifying the number of elements to send to each processor
      ! pmIssdis, integer, Scatterv Sendout DISpls
      !     * entry i specifies the displacement (relative to sendbuf rom which
      !     to take the outgoing data to process
      contains
! collects the size related info
      subroutine pmpi_master_gather_i(larr,n, nprocs, arr)
      integer, intent(in)                           :: n
      integer, dimension(n),intent(in)              :: larr
      integer, intent(in)                           :: nprocs
      integer, dimension(n*nprocs), intent(out)     :: arr

      integer, dimension(nprocs*n) :: temp
      !coarse grain communication
      call MPI_GATHER(larr, n, MPI_INTEGER, temp, n, MPI_INTEGER, MASTER,  MPI_COMM_WORLD, err)
      arr = RESHAPE( transpose( reshape(temp, (/n,nprocs/)) ), shape(arr) )
      end subroutine pmpi_master_gather_i

      
! broadcast the file name from master to all
! * file name
! * file name length
! * length without suffix
      subroutine pmpi_fname(infile, lth, trunc)
      integer, intent(inout)            :: lth
      character(len=lth), intent(inout) :: infile
      integer, intent(inout)            :: trunc

      call MPI_BCAST(lth, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, err)
      call MPI_BCAST(infile, lth, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, err)
      
      call MPI_BCAST(trunc, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, err)
      end subroutine pmpi_fname
      
      
      end module stmcomapi

