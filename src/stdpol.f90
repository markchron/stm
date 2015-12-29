      module stmdatpol
      use stmheader
	  integer 										:: nprocs, rank
      integer                                       :: nthread
	  integer 										:: stiflen, stiftrunc
	  character(BUF_LEN) 							:: stifnm ! input file name
	  character(BUF_LEN) 							:: stofnm ! input file name
	  integer, parameter 							:: FUNIT_OUT = 9
	  integer                                       :: FUNIT_LOG ! = 6?
	  character(BUF_LEN), dimension(10) 			:: titls_
	  
	  integer, dimension(0 : 100) 					:: stErrs
      ! integer control scalar variables 
	  integer, parameter 							:: NIcsv = 200
      integer, dimension(NIcsv), target         	:: stIcsv 
	  ! double control scalar variables 
      integer, parameter 							:: NDcsv = 100
      real(STDD), dimension(NDcsv), target      	:: stDcsv
      ! pointers
	  integer, pointer :: Nxd, Nyd, Nzd, Ngcll, Nwidx, Ngidx, Ngacl, & 
      Ncidx, Npidx, Neqn, Npeqn, Nxyplane, &
      Nlyid, Ngedges, Nevents, &
      fnstec, metric
	  integer, pointer :: Nlncs, Nlacs, Nlnbls, Nlnwls, Nlnlys, Nlneds
	  
	  ! stIgvs(szIdp), integer data pool; so far, only stIgvs(1:ofstIdp) was
	  ! assigned the pointers. 
	  integer 										:: szIdp, ofstIdp
      integer, dimension(:), allocatable, target    :: stIgvs
	  integer, dimension(:), pointer 				:: ictind, icsecd, &
	  ixly, iyly, izly, ialy, icidly, iwdir, ilyadj,  &
      ivadj, icdist, iwdist, iadjncy
	  ! ictind(Ngcll), integer, grid cell type index
	  ! 						normall type - all included in calculation
	  ! 						in-active 
	  ! 						pinch-out
      ! perforated layers info
      ! ixly(Nlyid), iyly(Nlyid), izly(Nlyid):UID of each perforated layer
      ! ialy(Nlyid) : connection of the layer: flow-from | flow-to
      ! icidly(Nlyid) : natural order of the perforated layer
      ! iwly(Nlyid) :
      ! iwdir(Nwidx) : well direction
      ! ilyadj(Nwidx+1) : perforated layers index of each well into perforated 
      ! layers info
      ! 
	  ! stDgvs(szDdp), double data pool; so far, only stDgvs(1:ofstDp) was
	  ! associated with pointers
	  integer 										:: szDdp, ofstDp
      real(STDD), dimension(:), allocatable, target :: stDgvs
	  real(STDD), dimension(:), pointer 			:: dcdx, dcdy, dcdz, ds3d, &
      dcareai, dcareaj, dcareak, area3d, dcpor, dcpermi, dcpermj, dcpermk, perm3d, &
      detrans, dethermtrs, & 
      dwllrad, dwllgeof, dwllfrac, dlylenfr, dlyskin, dlywitrs, dlywitherm, &
      dchtop
      ! dwllrad(Nwidx) : well radius
      ! dwllgeof(Nwidx) : geofactor
      ! dwllfrac(Nwidx) : well fraction 
      ! dlylenfr(Nlyid) : perforated length fraction
      ! dlyskin(Nlyid) : skin
      ! dlywitrs(Nlyid) : well index - flow
      ! dlywitherm    :               - heat conduction
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
	  Nlyid => stIcsv(13) 	! total perforated layers no.
	  Ngedges => stIcsv(14)  ! total connections (edges) no.
	  fnstec => stIcsv(15)  ! nine or five stencial
	  						! default, 5-points 
	  metric => stIcsv(16)  ! input units system
	  						! default, field units
      Nevents => stIcsv(17)  ! total events no.

	  Nlncs =>  stIcsv(7) 	! local grid cells number
	  Nlacs => 	stIcsv(8) 	! local active cells number
	  Nlnbls => stIcsv(9)   !
	  Nlnwls => stIcsv(4) 	! local well blks no.
	  Nlnlys => stIcsv(13)  ! local perforated layers no.
	  Nlneds => stIcsv(14)  ! local connections (edges) no.
	  
	  ofstIdp = 0; ofstDp = 0
	  fnstec = 0; metric = 0
	  end subroutine datpol_init_scalar
! after get the problem size (npt_init), initialize the array|vector pointers
	  subroutine datpol_init_vect
      integer           :: mxEdges
	  ! updates the connections number based on fnstec, Nxd, Nyd, Nzd	  
	  ! call numofconnects 
      mxEdges = estimate_no_connects(fnstec, Ngcll)
      ! properties on global grid cells			  
	  szIdp = Ngcll*3 + 5 * Nlyid + Ngcll+1 + Nwidx+1 + mxEdges + Nwidx

	  call palloc_i(stIgvs, szIdp, stErrs(1))
	  call setptr_i(stIgvs, szIdp, ofstIdp, ictind, Ngcll)
	  call setptr_i(stIgvs, szIdp, ofstIdp, icsecd, Ngcll)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ixly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, iyly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, izly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ialy, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, icidly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ivadj, Ngcll + 1)
	  call setptr_i(stIgvs, szIdp, ofstIdp, icdist, Ngcll)
	  call setptr_i(stIgvs, szIdp, ofstIdp, iwdist, Nwidx)

	  call setptr_i(stIgvs, szIdp, ofstIdp, iwdir, Nwidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ilyadj, Nwidx+1)
	  !call setptr_i(stIgvs, szIdp, ofstIdp, iadjncy, Ngedges)

	  szDdp = 10 * Ngcll + Nxyplane &
			  + 2 * mxEdges + 2 * Nlyid 
	  call palloc_d(stDgvs, szDdp, stErrs(2))
      call assoptr_d(stDgvs, szDdp, ofstDp, ds3d, 3*Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdx, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdy, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdz, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dchtop, Nxyplane)

      call assoptr_d(stDgvs, szDdp, ofstDp, area3d, 3*Ngcll)
      call setptr_d(stDgvs, szDdp, ofstDp, dcareai, Ngcll)
      call setptr_d(stDgvs, szDdp, ofstDp, dcareaj, Ngcll)
      call setptr_d(stDgvs, szDdp, ofstDp, dcareak, Ngcll)

	  call setptr_d(stDgvs, szDdp, ofstDp, dcpor, Ngcll)
      call assoptr_d(stDgvs, szDdp, ofstDp, perm3d, 3*Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcpermi, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcpermj, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcpermk, Ngcll)
      
      call setptr_d(stDgvs, szDdp, ofstDp, dwllrad, Nwidx)
      call setptr_d(stDgvs, szDdp, ofstDp, dwllgeof, Nwidx)
      call setptr_d(stDgvs, szDdp, ofstDp, dwllfrac, Nwidx)

      call setptr_d(stDgvs, szDdp, ofstDp, dlylenfr, Nlyid)
      call setptr_d(stDgvs, szDdp, ofstDp, dlyskin, Nlyid)
	  call setptr_d(stDgvs, szDdp, ofstDp, dlywitrs, Nlyid)
	  call setptr_d(stDgvs, szDdp, ofstDp, dlywitherm, Nlyid)

	  !call setptr_d(stDgvs, szDdp, ofstDp, detrans, Ngedges)
	  !call setptr_d(stDgvs, szDdp, ofstDp, dethermtrs, Ngedges)
	  end subroutine datpol_init_vect
! PURPOSE:
! based on the informations to updates the memory index
	  subroutine datpol_update_index
	  ! the active grid cells number
	  Ngacl = count( ictind > 0 )

      call set_geometry
      call set_wellbore
      call set_graph_connect
	  end subroutine datpol_update_index
! PURPOSE:
! calculate the geometry informations:
!      * area of each cell perpendicular three directions
      subroutine set_geometry
      use stmgeomet, only : set_area
      call set_area(dcdy, dcdz, dcareai)
      call set_area(dcdx, dcdz, dcareaj)
      call set_area(dcdx, dcdy, dcareak)
      end subroutine set_geometry
! PURPOSE:
! calculate the wellbore informations:
!   * uid (natural order) of perforated layer
      subroutine set_wellbore
      use stmwllbor, only : set_wll_cid, set_wlls_wi
      call set_wll_cid(Nlyid, ixly, iyly, izly, icidly, Nxd, Nyd) 
      call set_wlls_wi(Nwidx, Nlyid, ilyadj, icidly, iwdir, dwllrad, &
      dwllgeof, dwllfrac, dlyskin, dlylenfr, Ngcll, perm3d, ds3d,    &
      dlywitrs)
      end subroutine set_wellbore
! PURPOSE:
! estimate the connections number of Cartisian grid, based on 
! Nxd, Nyd, Nzd, fnstec
! ictind(Ngcll) [Grid cell type index with active grid cells info.] may reduce
! the connection lists
!      subroutine numofconnects
      ! accurately counts
      !use stmgeomet, only : no_totl_connections
      !call no_totl_connections( nxd, nyd, nzd, fnstec, Ngedges)
!      end subroutine numofconnects
      pure integer function estimate_no_connects(forn, nlls)
      integer, intent(in)           :: forn
      integer, intent(in)           :: nlls
      if(forn == 1 ) then
        estimate_no_connects = 11 * nlls
      else
        estimate_no_connects = 7 * nlls
      endif
      end function estimate_no_connects
!PURPOSE:
! update the graph connection info.
      subroutine set_graph_connect
      use stmgeomet, only : set_totl_connect_size, set_adjncy_trans
      call set_totl_connect_size(Nxd, Nyd, Nzd, fnstec, Ngedges, ivadj)
      call setptr_i(stIgvs, szIdp, ofstIdp, iadjncy, Ngedges)
      call setptr_d(stDgvs, szDdp, ofstDp, detrans, Ngedges)
      call set_adjncy_trans(Nxd,Nyd,Nzd,fnstec, Ngedges, ivadj, perm3d, area3d, ds3d, iadjncy, detrans)
      end subroutine set_graph_connect
! PURPOSE:
! domain partitioning algorithm
      subroutine set_partition_dist
      call set_dist_blocks
      end subroutine set_partition_dist
! PURPOSE:
! distribute the grid cells into each processor
      subroutine set_dist_clls
      use stmapimts, only : metis_api
      if (nprocs == 1) then
          icdist = 1
          iwdist = 1
          return
      endif
      ! give the graph of grid cells and weighted by transmissibility between grid cells
      call metis_api(0,rank,Ngcll, Ngedges, ivadj, iadjncy, detrans, nprocs, icdist)
      ! 1, metis_recursive_bisection 
      ! 0, metis_kway, keep contiguous of the partition 
      end subroutine set_dist_clls
! PURPOSE:
! distribute all the blocks (grid cell + wellblock) into each processor,
! the edge-cut weighted by transmissibility & well index.
!
! CSR-graph of grid (Ngcll, Ngedges, ivadj, iadjncy, detrans)
! CSR-graph of well (Nwidx, Nlyid, ilyadj, icidly, dlywitrs)
      subroutine set_dist_blocks
      use stmapimts, only : metis_api
      use stmgeomet, only : combine_cwll_adj
      use stmvtkxml, only : prt_csr
      integer           :: nvtxs, negds
      integer, dimension(Ngcll+Nwidx+1)         :: xadj
      integer, dimension(Ngedges+2*Nlyid)       :: adjncy
      real(STDD), dimension(Ngedges+2*Nlyid)    :: adjwgt
      integer, dimension(Ngcll+Nwidx)           :: xdist

      nvtxs = Ngcll + Nwidx
      negds = Ngedges + 2 * Nlyid

      call combine_cwll_adj(Ngcll, Ngedges, ivadj, iadjncy, detrans, &
      Nwidx, Nlyid, ilyadj, icidly, dlywitrs,  &
      nvtxs, negds, xadj, adjncy, adjwgt)
      
      call prt_csr(FUNIT_OUT, nvtxs, negds, xadj, adjncy, adjwgt,"trans-wellindex", 1)
      
      if(nprocs == 1) then
          xdist = 1
      else
        call metis_api(0,rank,nvtxs, negds, xadj, adjncy, adjwgt, nprocs, xdist)
        ! 1, metis_recursive_bisection 
        ! 0, metis_kway, failed !!!
      endif

      icdist = xdist(1:Ngcll)
      iwdist = xdist(Ngcll + 1 : Ngcll + Nwidx)
      end subroutine set_dist_blocks

! PURPOSE:
! de-associate the pointers
! release the memory
      subroutine datpol_free
      NULLIFY(ictind, icsecd)
	  NULLIFY(ixly, iyly, izly, ialy, icidly, iwdir, ilyadj, ivadj, iadjncy)
	  NULLIFY(dcdx, dcdy, dcdz, dchtop, dcpor, dcpermi, dcpermj, dcpermk)
	  NULLIFY(detrans, dethermtrs)
      NULLIFY(dwllrad, dwllgeof, dwllfrac)
      NULLIFY(dlylenfr, dlyskin, dlywitrs, dlywitherm)
	  deallocate(stDgvs, STAT = stErrs(11))
	  deallocate(stIgvs, STAT = stErrs(12))
	  end subroutine datpol_free
! count the maximum columns no. for CSR matrix
      pure integer function maxcol_csr(n, ai)
      integer, intent(in)           :: n
      integer, dimension(n+1), intent(in) :: ai
      integer :: i, ncol
      maxcol_csr = ai(2) - ai(1) 
      do i= 2, n
        ncol = ai(i+1) - ai(i)
        if(maxcol_csr < ncol) maxcol_csr = ncol
      enddo
      end function maxcol_csr
      end module stmdatpol
