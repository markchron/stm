	  module stmdatnpt
      use stmheader
	  use stmdatpol

	  contains
! pre-read the data size | structure related info.
! predefine all these  values by revise the program
      subroutine npt_init
      Nxd = 9
	  Nyd = 9
      Nzd = 4
	  Nwidx = 3
	  Nlyidx = 3 * Nzd  ! defined after go through the datset, here, 3 vertical wells
	  nfstec = 0 		! pre-read:  discretized format, 5 or 9 points
	  metric = 0 		! field units

	  Ncidx = 2
	  Npidx = 3
	  call npt_update
	  call datpol_init_vect
      end subroutine npt_init
! update the data size related info 
	  subroutine npt_update
	  Nxyplane = Nxd * Nyd
      Ngcll = Nxyplane * Nzd
	  Ngidx = Ngcll + Nwidx

	  Neqn = Ncidx + 1 + Npidx + 1
	  Npeqn = Ncidx + 1
	  end subroutine npt_update

	  subroutine npt_set_properties
      integer :: i, k
      ictind = (/ ( 5, 2, 2, 2, 2, 2, 2, 2, 5,  &
			  0, 4, 1, 1, 1, 1, 1, 3, 0, &
			  0, 0, 4, 1, 1, 1, 3, 0, 0,  &
			  0, 0, 0, 4, 1, 3, 0, 0, 0, &
			  0, 0, 0, 0, 6, 0, 0, 0, 0, &
			  (0, i=1, 36), k=1, Nzd) /)
	  dcdx = 29.17 ! ft
	  dcdy = 29.17 ! ft
	  dcdz = (/ (10., i=1, Nxyplane),  (20., i=1, Nxyplane), &
			  (25., i=1, Nxyplane),  (25., i=1, Nxyplane)  /) ! ft, KVAR
	  dchtop = 1500 ! ft
	  
	  dcpor = 0.3
	  dcpermi  =  (/ (2000., i=1, Nxyplane),  (500., i=1, Nxyplane), &
			  (1000., i=1, Nxyplane),  (2000., i=1, Nxyplane)  /) ! Darcy, KVAR
	  dcpermj = dcpermi
	  dcpermk = dcpermi * 0.5
	  end subroutine npt_set_properties
	  end module stmdatnpt
