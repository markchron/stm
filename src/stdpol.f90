      module stmdatapool
      use stmheader
      integer, parameter :: stNIntege = 100
      integer, parameter :: stNDouble = 100
      integer, dimension(stNIntege), target         :: stIgss
      integer, dimension(:), allocatable, target    :: stIgvs
      real(STDD), dimension(stNDouble), target      :: stDgss
      real(STDD), dimension(:), allocatable, target :: stDgvs

      ! pointers
      integer, pointer :: nxd, nyd, nzd, ncidx, npidx

      end module stmdatapool
