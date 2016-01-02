      module stmheader
      implicit none
      integer, parameter :: STDD 			= 8
      integer, parameter :: BUF_LEN			= 72
      real(kind=STDD), parameter :: ST_PI = 3.14159265359
      real(kind=STDD), parameter :: ST_TOL_EPSILON = 1.e-13
      real(kind=STDD), parameter :: ST_FPE_TOL = 1.e-13
      real(kind=STDD), parameter :: ST_FPE_INF = 9.e30
      
      integer, parameter         :: MASTER = 0
      integer, parameter         :: ST_OMP_LIMIT = 5000
      
       ! associate array pointer to a target data pool
      interface set_ptr
        module procedure set_ptr_d
        module procedure set_ptr_i
      end interface  ! set_ptr   
      contains
! PURPOSE:
! dynamic memory, check the allocation 
      subroutine palloc_i(arr, n, err)
      implicit none
      integer, dimension(:),allocatable, intent(inout) :: arr
      integer, intent(in) :: n
      integer, intent(out) :: err
      if(allocated(arr)) then
        err = 1
        deallocate(arr)
      endif
      allocate(arr(n),STAT=err)
      end subroutine palloc_i
      subroutine palloc_d(arr, n, err)
      implicit none
      real(STDD), dimension(:),allocatable, intent(inout) :: arr
      integer, intent(in) :: n
      integer, intent(out) :: err
      if(allocated(arr)) then
        err = 1
        deallocate(arr)
      endif
      allocate(arr(n),STAT=err)
      end subroutine palloc_d
     
      subroutine setptr_i(pool, n, off, ptr, np)
      implicit none
      integer, intent(in) :: n
      integer, dimension(n), target, intent(inout) :: pool
      integer, intent(inout) :: off
      integer, intent(in) :: np
      integer, dimension(:), pointer :: ptr
      call assoptr_i(pool,n,off,ptr,np) 
      off  = off + np
      end subroutine setptr_i
      subroutine assoptr_i(pool, n, off, ptr, np)
      implicit none
      integer, intent(in) :: n
      integer, dimension(n), target, intent(inout) :: pool
      integer, intent(in) :: off
      integer, intent(in) :: np
      integer, dimension(:), pointer :: ptr
      integer       :: offB 
      integer       :: offA
      offA = off + 1
      offB = off + np
      ptr  => pool(offA : offB : 1)
      end subroutine assoptr_i
      subroutine setptr_d(pool, n, off, ptr, np)
      implicit none
      integer, intent(in) :: n
      real(STDD), dimension(n), target, intent(inout) :: pool
      integer, intent(inout) :: off
      integer, intent(in) :: np
      real(STDD), dimension(:), pointer :: ptr
      call assoptr_d(pool, n, off, ptr, np)
      off  = off + np
      end subroutine setptr_d
      subroutine assoptr_d(pool, n, off, ptr, np)
      implicit none
      integer, intent(in) :: n
      real(STDD), dimension(n), target, intent(inout) :: pool
      integer, intent(in) :: off
      integer, intent(in) :: np
      real(STDD), dimension(:), pointer :: ptr
      integer       :: offB 
      integer       :: offA
    
      offA = off + 1
      offB = off + np
      ptr  => pool(offA : offB : 1)
      end subroutine assoptr_d

      subroutine set_ptr_d(pool, n, off, ptr, np, RESET)
      integer, intent(in) :: n
      real(kind=STDD), dimension(n), target, intent(inout) :: pool
      integer, intent(out) :: off
      integer, intent(in) :: np
      real(kind=STDD), dimension(:), pointer :: ptr
      logical, optional :: RESET
      integer, save :: offB = 0
      integer       :: offA

      if(PRESENT(RESET).and.RESET) offB = 0
      
      off  = offB
      offA = offB + 1
      offB = offB + np
      ptr  => pool(offA : offB : 1)      
      end subroutine set_ptr_d
      subroutine set_ptr_i(pool, n, off, ptr, np, RESET)
      implicit none
      integer, intent(in) :: n
      integer, dimension(n), target, intent(inout) :: pool
      integer, intent(out) :: off
      integer, intent(in) :: np
      integer, dimension(:), pointer :: ptr
      logical, optional :: RESET
      integer, save :: offB = 0
      integer       :: offA

      if(PRESENT(RESET).and.RESET) offB = 0
      
      off  = offB
      offA = offB + 1
      offB = offB + np
      ptr  => pool(offA : offB : 1)
      end subroutine set_ptr_i
! PURPOSE:
! convert CSR into CSC format B=A, A in CSR and B in CSC matrix
! 
! input:
! integer nr          -- number of rows in 2D density format
! integer nc          -- number of columns in density format
! integer ap(nr+1)    -- row  pointer in CSR
! integer aj(nnz(A))  -- column indices
! real av(nnz(A))     -- non zero values
!
! output:
! integer bp(nc+1)    -- column pointer
! integer bi(nnz(A))  -- row indices
! real bv(nnz(A))     -- non zero values
!
! note
! bp, bi, bx must be preallocated
      subroutine csr2csc(nr, nc,nnz, ap, aj, av, bp, bi, bv)
      implicit none
      integer, intent(in)                       :: nr
      integer, intent(in)                       :: nnz
      integer, dimension(nr+1), intent(in)                  :: ap
      integer, dimension(nnz), intent(in)                   :: aj
      real(kind=STDD), dimension(nnz), intent(in), optional :: av
      integer, intent(in)                       :: nc
      integer, dimension(nc+1), intent(out)                 :: bp
      integer, dimension(nnz), intent(out)                  :: bi
      real(kind=STDD), dimension(nnz), intent(out), optional:: bv
      integer :: cumsum, temp
      integer :: l, r, c, dest
      ! compute number of non-zero entries per column of density matrix
      bp = 0
!$OMP PARALLEL DO PRIVATE(l) reduction(+:bp)
      do l=1, nnz
        bp(aj(l)) = bp(aj(l)) + 1
      enddo
!$OMP END PARALLEL DO

      ! cummulate sum the number of non-zero entries per column to get bp()
      cumsum = 1
      ! C/C++ , cumsum starts with 0, Fortran starts with 1
      do c=1, nc
        temp = bp(c)
        bp(c) = cumsum
        cumsum = cumsum + temp
      enddo
      bp(nc+1) = cumsum  ! == nnz + 1
                         ! bi(bp(c) : bp(c+1)-1) is the column

      ! go through csr structure once more, fill in output matrix
      do r=1, nr
        do l=ap(r), ap(r+1) - 1
          c = aj(l)
          dest = bp(c)

          bi(dest) = r
          if(PRESENT(av).and.PRESENT(bv)) then
            bv(dest) = av(l)
          endif
          ! shift forward bp(:) 
          bp(c) = bp(c) + 1
        enddo ! l
      enddo! r

      ! shift back bp(:)
      do c= nc+1, 2, -1
        bp(c) = bp(c-1)
      enddo
      bp(1) = 1
      end subroutine csr2csc
! add up the array with the previous entries, each entry has same weigth
      subroutine nums2displs_i(n, nums, displs)
      integer, intent(in)                   :: n
      integer, dimension(n), intent(in)     :: nums
      integer, dimension(n), intent(out)     :: displs
      integer    :: l
      displs(1) = 0
      forall (l = 1 : n-1)
      displs(l+1) = displs(l) + nums(l)
      end forall
      end subroutine nums2displs_i
      subroutine nums2displs_d(n, nums, displs)
      integer, intent(in)                      :: n
      real(STDD), dimension(n), intent(in)     :: nums
      real(STDD), dimension(n), intent(out)     :: displs
      integer    :: l
      displs(1) = 0.d0
      forall (l = 1 : n-1)
      displs(l+1) = displs(l) + nums(l)
      end forall
      end subroutine nums2displs_d

      end module stmheader
