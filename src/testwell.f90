      program testwll
      use stmheader
      use stmwllbor
      integer :: n
      real(STDD), dimension(:), allocatable, target :: memo
      real(STDD), pointer :: wfrac, rw, skin, geofac
      real(STDD), dimension(:), pointer :: permi, permj, dx, dy, dz, hf, &
      wi

      n = 4
      call intial_memo
      call set_memo
      call wll_wi(wfrac, geofac, rw, skin, n, permi, permj, dz, dx, dy, hf, wi)

      contains
      subroutine intial_memo
      integer :: err, nmemo, of

      nmemo = 7 * n + 4
      call palloc_d(memo, nmemo, err)
      wfrac     => memo(1) 
      rw        => memo(2) 
      skin      => memo(3) 
      geofac    => memo(4)

      of = 4
      
      call setptr_d(memo, nmemo, of, permi, n)
      call setptr_d(memo, nmemo, of, permj, n)
      call setptr_d(memo, nmemo, of, dx, n)
      call setptr_d(memo, nmemo, of, dy, n)
      call setptr_d(memo, nmemo, of, dz, n)
      call setptr_d(memo, nmemo, of, hf, n)
      call setptr_d(memo, nmemo, of, wi, n)
      end subroutine intial_memo

      subroutine set_memo
      wfrac = 0.25
      rw = 0.3 ! ft
      skin = 0 
      geofac = 0.249

      permi = (/ 2000, 500, 1000, 2000 /) ! mD
      permj = permi
!      dx  = 29.17       ! ft
!      dy  = 29.17       ! ft
      dx  = 15.44       ! ft
      dy  = 15.44       ! ft
      dz  = (/10 , 20, 25, 25 /) ! ft
      hf = 1  ! well perforated length fraction
      end subroutine set_memo
      end program testwll
