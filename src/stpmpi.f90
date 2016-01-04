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
      subroutine pmpi_master_gatherv_i(ns, sendbuf, sendcount,          &
      nr, recvbuf, np, recvcounts, recvdisp)
      integer, intent(in)                           :: ns
      integer, dimension(ns), intent(in)            :: sendbuf
      integer, intent(in)                           :: sendcount
      integer, intent(in)                           :: nr
      integer, dimension(nr), intent(out)           :: recvbuf
      integer, intent(in)                           :: np
      integer, dimension(np), intent(in)            :: recvcounts
      integer, dimension(np), intent(in)            :: recvdisp

      call MPI_GATHERV(sendbuf, sendcount, MPI_INTEGER, recvbuf,        &
      recvcounts, recvdisp, MPI_INTEGER, MASTER, MPI_COMM_WORLD, err)
      end subroutine pmpi_master_gatherv_i
! each process sends distinct data to each of the receives. The j-th
! block sent from process i is received by process j and is placed in
! the i-th block of recvbuf
! The type of data must be consisitent. The amount of data sent must be
! equal to the amount of data received, pairwise between every pair of
! processes. 
      subroutine pmpi_all_exchange_i(n, sendbuf, sendcount, recvbuf, recvcount)
      integer, intent(in)                           :: n
      integer, dimension(n), intent(in)             :: sendbuf
      integer, intent(in)                           :: sendcount
      integer, dimension(n), intent(out)            :: recvbuf
      integer, intent(in)                           :: recvcount
      
      call MPI_ALLTOALL(sendbuf, sendcount, MPI_INTEGER, recvbuf,       &
      recvcount, MPI_INTEGER, MPI_COMM_WORLD, err)
      end subroutine pmpi_all_exchange_i
      subroutine pmpi_all_exchangev_i(ns, sendbuf, np, sendcounts, senddisp,&
      nr, recvbuf, recvcounts, recvdisp)
      integer, intent(in)                           :: ns
      integer, dimension(ns), intent(in)            :: sendbuf
      integer, intent(in)                           :: np
      integer, dimension(np), intent(in)            :: sendcounts
      integer, dimension(np), intent(in)            :: senddisp
      integer, intent(in)                           :: nr
      integer, dimension(nr), intent(out)            :: recvbuf
      integer, dimension(np), intent(in)            :: recvcounts
      integer, dimension(np), intent(in)            :: recvdisp
      
      call MPI_ALLTOALLV(sendbuf, sendcounts, senddisp, MPI_INTEGER,    &
      recvbuf, recvcounts, recvdisp, MPI_INTEGER, MPI_COMM_WORLD, err)
      end subroutine pmpi_all_exchangev_i
      
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

