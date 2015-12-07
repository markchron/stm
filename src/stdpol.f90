      module stmdatpol
      use stmheader
	  character(BUF_LEN) 							:: stifnm
	  character(BUF_LEN), dimension(10) 			:: titls_
	  integer, dimension(100) 						:: stErrs
      ! integer control scalar variables 
	  integer, parameter 							:: NIcsv = 200
      integer, dimension(NIcsv), target         	:: stIcsv 
	  ! double control scalar variables 
      integer, parameter 							:: NDcsv = 100
      real(STDD), dimension(NDcsv), target      	:: stDcsv
      ! pointers
	  integer, pointer :: Nxd, Nyd, Nzd, Ngcll, Nwidx, Ngidx, Ngacl, Ncidx, Npidx, Neqn, Npeqn, Nxyplane, &
			  Nlyidx, fnstec, metric, Ngedges
	  integer, pointer :: Nlncs, Nlacs, Nlnbls, Nlnwls, Nlnlys, Nlneds
	  
	  ! stIgvs(szIdp), integer data pool; so far, only stIgvs(1:ofstIdp) was
	  ! assigned the pointers. 
	  integer 										:: szIdp, ofstIdp
      integer, dimension(:), allocatable, target    :: stIgvs
	  integer, dimension(:), pointer 				:: ictind, icsecd, &
	  ixly, iyly, izly, ialy, iwly, ivadj, iadjncy
	  ! ictind(Ngcll), integer, grid cell type index
	  ! 						normall type - all included in calculation
	  ! 						in-active 
	  ! 						pinch-out
	  ! stDgvs(szDdp), double data pool; so far, only stDgvs(1:ofstDp) was
	  ! associated with pointers
	  integer 										:: szDdp, ofstDp
      real(STDD), dimension(:), allocatable, target :: stDgvs
	  real(STDD), dimension(:), pointer 			:: dcdx, dcdy, dcdz, dcpor, &
	  dchtop, dcpermi, dcpermj, dcpermk, &
      detrans, dethermtrs, dwlltrs, dwllthermtrs
	  contains
! initialize the scalar pointers to the global datapool.
	  subroutine datpol_init_scalar
      Nxd => 	stIcsv(1) 	! Nxd, N in x
	  Nyd => 	stIcsv(2) 	! Nyd, N in y
	  Nzd => 	stIcsv(3) 	! Nzd, N in z
	  Nwidx => 	stIcsv(4)  	! Nwidx, total wells number
	  Ncidx => 	stIcsv(5) 	! Ncidx, total components no.
	  Npidx => 	stIcsv(6) 	! Npidx, total phases no.

	  Ngcll => 	stIcsv(7) 	! total grid cells number
	  Ngacl => 	stIcsv(8) 	! total active grid cells number
      Ngidx =>  stIcsv(9) 	! total global blocks number
	  Npeqn => 	stIcsv(10)  ! primary eqns no.
	  Neqn 	=> 	stIcsv(11)  ! Ncidx + 1 + Npidx + 1 + Nequil
	  Nxyplane => stIcsv(12) ! Nxd * Nyd
	  Nlyidx => stIcsv(13) 	! total perforated layers no.
	  Ngedges => stIcsv(14)  ! total connections (edges) no.
	  fnstec => stIcsv(15)  ! nine or five stencial
	  						! default, 5-points 
	  metric => stIcsv(16)  ! input units system
	  						! default, field units

	  Nlncs =>  stIcsv(7) 	! local grid cells number
	  Nlacs => 	stIcsv(8) 	! local active cells number
	  Nlnbls => stIcsv(9)   !
	  Nlnwls => stIcsv(4) 	! local well blks no.
	  Nlnlys => stIcsv(13)  ! local perforated layers no.
	  Nlneds => stIcsv(14)  ! local connections (edges) no.
	  
	  ofstIdp = 0; ofstDp = 0
	  nfstec = 0; metric = 0
	  end subroutine datpol_init_scalar
! after get the problem size (npt_init), initialize the array|vector pointers
	  subroutine datpol_init_vect
	  ! updates the connections number based on nfstec, Nxd, Nyd, Nzd	  
	  call numofconnects 
      ! properties on local grid cells			  
	  szIdp = Ngcll *  2 + 5 * Nlyidx + Ngcll+1 + Ngedges

	  call palloc_i(stIgvs, szIdp, stErrs(1))
	  call setptr_i(stIgvs, szIdp, ofstIdp, ictind, Ngcll)
	  call setptr_i(stIgvs, szIdp, ofstIdp, icsecd, Ngcll)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ixly, Nlyidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, iyly, Nlyidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, izly, Nlyidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ialy, Nlyidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, iwly, Nlyidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ivadj, Ngcll + 1)
	  call setptr_i(stIgvs, szIdp, ofstIdp, iadjncy, Ngedges)

	  szDdp = 7 * Ngcll + Nxyplane &
			  + 2 * Ngedges + 2 * Nlyidx 
	  call palloc_d(stDgvs, szDdp, stErrs(2))
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdx, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdy, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdz, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dchtop, Nxyplane)

	  call setptr_d(stDgvs, szDdp, ofstDp, dcpor, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcpermi, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcpermj, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcpermk, Ngcll)

	  call setptr_d(stDgvs, szDdp, ofstDp, detrans, Ngedges)
	  call setptr_d(stDgvs, szDdp, ofstDp, dethermtrs, Ngedges)
	  call setptr_d(stDgvs, szDdp, ofstDp, dwlltrs, Nlyidx)
	  call setptr_d(stDgvs, szDdp, ofstDp, dwllthermtrs, Nlyidx)

	  end subroutine datpol_init_vect

	  subroutine datpol_update_index
      integer :: i, n
	  n = 0
	  do i=1, Nlncs
        if(ictind(i) > 0 ) n = n + 1
      enddo

	  n = count( ictind > 0 )
	  Nlacs = n
	  end subroutine datpol_update_index
      subroutine numofconnects
      Ngedges = 5 * Ngcll
      end subroutine numofconnects
      end module stmdatpol
