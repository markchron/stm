      module stmheader
      implicit none
      integer, parameter :: STDD 			= 8
      integer, parameter :: BUF_LEN			= 72
      real(kind=STDD), parameter :: ST_PI = 3.14159265359
      real(kind=STDD), parameter :: ST_TOL_EPSILON = 1.e-13
      real(kind=STDD), parameter :: ST_FPE_TOL = 1.e-13
      real(kind=STDD), parameter :: ST_FPE_INF = 9.e30
      
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
      integer       :: offB 
      integer       :: offA
    
      offA = off + 1
      offB = off + np
      ptr  => pool(offA : offB : 1)
      off  = offB
      end subroutine setptr_i
      subroutine setptr_d(pool, n, off, ptr, np)
      implicit none
      integer, intent(in) :: n
      real(STDD), dimension(n), target, intent(inout) :: pool
      integer, intent(inout) :: off
      integer, intent(in) :: np
      real(STDD), dimension(:), pointer :: ptr
      integer       :: offB 
      integer       :: offA
    
      offA = off + 1
      offB = off + np
      ptr  => pool(offA : offB : 1)
      off  = offB
      end subroutine setptr_d      

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
      end module stmheader
