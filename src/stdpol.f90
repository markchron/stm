      module stmdatpol
      use stmheader
      integer                                       :: nprocs, rank, deck
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
      ! local 
      integer, pointer :: Nlncs, Nlacs, Nlnbls, Nlnwls, Nlnlys, Nlneds, &
      Nlbinn, Nlbbrd, Nlbext, Nlbint, Nlcinn, Nlcbrd, Nlcext, Nlcint,   &
      Nlwinn, Nlwbrd, Nlwext, Nlwint
      ! Nlnbls, local blocks (cells + wells) no. (Nlncs+Nlnwls)
      !                             = Nlbinn + Nlbbrd + Nlbext
      !                             = Nlbint + Nlbext
      ! Nlncs, local grid cells no. (internal+external cells)
      !                             = Nlcint + Nlcext
      ! Nlnwls, local well blks no. (internal + external wells)
      !                             = Nlwint + Nlwext
      ! Nlacs, local active cells no. 
      ! Nlnlys, local perforated layers no.
      ! Nlneds, local connections no.

      ! Nlbinn, local inner blocks no.
      ! Nlbbrd, local border blocks no.
      ! Nlbext, local external blocks no.
      ! Nlbint, internal blocks no. = Nlbinn + Nlbbrd

	  
	  ! stIgvs(szIdp), integer data pool; so far, only stIgvs(1:ofstIdp) was
	  ! assigned the pointers. 
	  integer 										:: szIdp, ofstIdp
      integer, dimension(:), allocatable, target    :: stIgvs
      integer, dimension(:), pointer                :: ictind, icsecd, icgcat,&
      icdist, iccatl, iclcid, &
      iwdist, iwdir, iwcatl, iwlcid, &
      ixly, iyly, izly, ialy, icidly, &
      ivadj, ilyadj, iadjncy, &
      ibdist, ibcatl, iblcid
	  ! ictind(Ngcll), integer, grid cell type index
	  ! 						normall type - all included in calculation
	  ! 						in-active 
	  ! 						pinch-out
      ! icsecd(Ngcll)
      ! iccatl(Ngcll), integer, grid cell category on local partition
      !                         = 0, NOT ON current deck
      !                         = 1, inner
      !                         = 2, border
      !                         = 3, external
      ! icgcat(Ngcll), integer, grid cell location info.
      !                         BTEST(icgcat, 1) I-
      !                         BTEST(icgcat, 2) I+
      !                         BTEST(icgcat, 3) J-
      !                         BTEST(icgcat, 4) J+
      !                         BTEST(icgcat, 5) underburden, kdir up:K-,down:K+
      !                         BTEST(icgcat, 6) overburden, kdir up:K+,down:K-
      ! icdist(Ngcll), integer, grid cell assigned domain|deck index
      ! 
      ! iwdist(Nwidx), integer, wellbore block assigned domain|deck index
      ! iwdir(Nwidx) : well direction
      ! perforated layers info
      ! ixly(Nlyid), iyly(Nlyid), izly(Nlyid):UID of each perforated layer
      ! ialy(Nlyid) : connection of the layer: flow-from | flow-to
      ! icidly(Nlyid) : natural order of the perforated layer
      ! ilyadj(Nwidx+1) : perforated layers index of each well into perforated 
      ! layers info
      ! 
      ! iblcid(Ngidx), integer, local index (inner/border/external) of all the
      ! blocks on current subdeck
      ! iclcid(Ngcll),  cells
      ! iwlcid(Nwidx),  wells
      ! 
	  ! stDgvs(szDdp), double data pool; so far, only stDgvs(1:ofstDp) was
	  ! associated with pointers
	  integer 										:: szDdp, ofstDp
      real(STDD), dimension(:), allocatable, target :: stDgvs
	  real(STDD), dimension(:), pointer             :: dcdx, dcdy, dcdz, ds3d, &
      dcareai, dcareaj, dcareak, area3d, dcpor, dcpermi, dcpermj, dcpermk, perm3d, &
      detrans, dethermtrs, & 
      dwllrad, dwllgeof, dwllfrac, &
      dlylenfr, dlyskin, dlywitrs, dlywitherm, &
      dchtop, devnts
      ! dwllrad(Nwidx) : well radius
      ! dwllgeof(Nwidx) : geofactor
      ! dwllfrac(Nwidx) : well fraction 
      ! dlylenfr(Nlyid) : perforated length fraction
      ! dlyskin(Nlyid) : skin
      ! dlywitrs(Nlyid) : well index - flow
      ! dlywitherm    :               - heat conduction

      ! devnts(Events), date|events info.
      !--------------------------------------------------------------------
      ! mpi datpol
      !--------------------------------------------------------------------
      integer                                    :: pmnIcs
      integer                                    :: pmofIc
      integer, dimension(:), allocatable, target :: pmIcsv
      ! pmIcsv(:), Parallel Message passing Integer Control System Vector
      integer, dimension(:), pointer            ::  pmInec, pmIdisec, &
      pmInew, pmIdisew,  &
      pmInmas, pmInupcl, pmInupwl, pmInupbk, pmInlcll, pmInlwll, pmInlblk
      ! pmInmas(6*nprocs), includes the following parts
      ! pmInupcl(nprocs)        :: internal cells
      ! pmInupwl(nprocs)        :: internal wells
      ! pmInupbk(nprocs)        :: internal blocks
      ! pmInlcll(nprocs)        :: local cells no. of each subdeck
      ! pmInlwll(nprocs)        :: local wells no. of each subdeck
      ! pmInlblk(nprocs)        :: local blocks no. 
      !--------------------------------------------------------------------
      ! local datapool
      !--------------------------------------------------------------------
      integer                                       :: szIloc, ofIloc
      integer, dimension(:), allocatable, target    :: dpIloc
      integer, dimension(:), pointer                :: ibgid, icgid, iwgid
      ! ibgid, block natural order
      ! icgid, cell natural order
      ! iwgid, well natural order
      
      
      contains
! set up the all2all index
      subroutine set_exchange_index
      use stmgeomet, only : set_locid, set_numdisp_extofadjdeck
      ! order the local vertex
      ! *cell
      call set_locid(Ngcll, iccatl, icdist, Nlcinn, Nlcbrd, Nlcext, & 
      iclcid, Nlncs, icgid) 
      call set_numdisp_extofadjdeck(Ngcll, icdist, iccatl, nprocs, pmInec, pmIdisec)
      ! *well
      call set_locid(Nwidx, iwcatl, iwdist, Nlwinn, Nlwbrd, Nlwext, &
      iwlcid, Nlnwls, iwgid)
      call set_numdisp_extofadjdeck(Nwidx, iwdist, iwcatl, nprocs, pmInew, pmIdisew)

      end subroutine set_exchange_index
! initialize the scalar pointers to the global datapool.
	  subroutine datpol_init_scalar
      Nxd => 	stIcsv(1) 	! Nxd, N in x
	  Nyd => 	stIcsv(2) 	! Nyd, N in y
	  Nzd => 	stIcsv(3) 	! Nzd, N in z
	  Ngacl => 	stIcsv(4) 	! total active grid cells number
	  Ncidx => 	stIcsv(5) 	! Ncidx, total components no.
	  Npidx => 	stIcsv(6) 	! Npidx, total phases no.

	  Ngcll => 	stIcsv(7) 	! total grid cells number
	  Nwidx => 	stIcsv(8)  	! Nwidx, total wells number
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

      ofstIdp = 0; ofstDp = 0; pmofIc = 0; ofIloc = 0
      fnstec = 0; metric = 0
      end subroutine datpol_init_scalar
! PURPOSE:
      subroutine datpol_init_scalar_loc
      Nlcinn    => stIcsv(101) ! cells
      Nlcbrd    => stIcsv(104)
      Nlcext    => stIcsv(107)
      Nlcint    => stIcsv(110)
      Nlncs     => stIcsv(113) ! local cells

      Nlwinn    => stIcsv(102) ! wells
      Nlwbrd    => stIcsv(105)
      Nlwext    => stIcsv(108)
      Nlwint    => stIcsv(111)
      Nlnwls    => stIcsv(114) ! local well 

      Nlbinn    => stIcsv(103) ! inner blks
      Nlbbrd    => stIcsv(106) ! border blks
      Nlbext    => stIcsv(109) ! external blks
      Nlbint    => stIcsv(112) ! internal blks.
      Nlnbls    => stIcsv(115) ! local blks. = internal + external

      Nlacs     => stIcsv(116)
      Nlnlys    => stIcsv(117)
      Nlneds    => stIcsv(118)

      end subroutine datpol_init_scalar_loc
! after get the problem size (npt_init). get memory allocation
      subroutine datpol_init_memo
      call datpol_init_vect
      call pmpi_init_vect
      end subroutine datpol_init_memo
      
! initialize the array|vector pointers
	  subroutine datpol_init_vect
      integer           :: mxEdges
	  ! updates the connections number based on fnstec, Nxd, Nyd, Nzd	  
	  ! call numofconnects 
      mxEdges = estimate_no_connects(fnstec, Ngcll)
      ! properties on global grid cells			  
	  szIdp = Ngcll*6           &
      + Nwidx * 4               &
      + Nlyid * 5               &
      + Ngcll+1 + Nwidx+1       &
      + mxEdges 

      call palloc_i(stIgvs, szIdp, stErrs(1))
      call setptr_i(stIgvs, szIdp, ofstIdp, icsecd, Ngcll)
      call setptr_i(stIgvs, szIdp, ofstIdp, ictind, Ngcll)
      call setptr_i(stIgvs, szIdp, ofstIdp, icgcat, Ngcll)

      call assoptr_i(stIgvs, szIdp, ofstIdp, ibdist, Ngidx)
      call setptr_i(stIgvs, szIdp, ofstIdp, icdist, Ngcll)
      call setptr_i(stIgvs, szIdp, ofstIdp, iwdist, Nwidx)

      call assoptr_i(stIgvs, szIdp, ofstIdp, ibcatl, Ngidx)
      call setptr_i(stIgvs, szIdp, ofstIdp, iccatl, Ngcll)
      call setptr_i(stIgvs, szIdp, ofstIdp, iwcatl, Nwidx)

	  call setptr_i(stIgvs, szIdp, ofstIdp, iwdir, Nwidx)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ixly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, iyly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, izly, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ialy, Nlyid)
	  call setptr_i(stIgvs, szIdp, ofstIdp, icidly, Nlyid)

	  call setptr_i(stIgvs, szIdp, ofstIdp, ivadj, Ngcll + 1)
	  call setptr_i(stIgvs, szIdp, ofstIdp, ilyadj, Nwidx+1)
	  !call setptr_i(stIgvs, szIdp, ofstIdp, iadjncy, Ngedges)

      call assoptr_i(stIgvs, szIdp, ofstIdp, iblcid, Ngidx)
      call setptr_i(stIgvs, szIdp, ofstIdp, iclcid, Ngcll)
      call setptr_i(stIgvs, szIdp, ofstIdp, iwlcid, Nwidx)

      szDdp = Ngcll * 10            &
      + Nwidx * 3                   &
      + Nlyid * 4                   &
      + Nxyplane + Nevents          &
      + mxEdges * 2                 

	  call palloc_d(stDgvs, szDdp, stErrs(2))

      call assoptr_d(stDgvs, szDdp, ofstDp, ds3d, 3*Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdx, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdy, Ngcll)
	  call setptr_d(stDgvs, szDdp, ofstDp, dcdz, Ngcll)
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
	  call setptr_d(stDgvs, szDdp, ofstDp, dchtop, Nxyplane)
      call setptr_d(stDgvs, szDdp, ofstDp, devnts, Nevents)
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
!      * grid cell location, BITWISE info. 
      subroutine set_geometry
      use stmgeomet, only : set_area, set_clls_gcate
      ! area
      call set_area(dcdy, dcdz, dcareai)
      call set_area(dcdx, dcdz, dcareaj)
      call set_area(dcdx, dcdy, dcareak)
      ! location info.
      call set_clls_gcate(Nxd, Nyd, Nzd, icgcat)
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
      ! counts total edges no. based on grid cells and five|nine stencial
      call set_totl_connect_size(Nxd, Nyd, Nzd, fnstec, Ngedges, ivadj)
      ! dynamic allocation
      call setptr_i(stIgvs, szIdp, ofstIdp, iadjncy, Ngedges)
      call setptr_d(stDgvs, szDdp, ofstDp, detrans, Ngedges)
      call setptr_d(stDgvs, szDdp, ofstDp, dethermtrs, Ngedges)

      call set_adjncy_trans(Nxd,Nyd,Nzd,fnstec, Ngedges, ivadj, perm3d, area3d, ds3d, iadjncy, detrans)
      end subroutine set_graph_connect
! PURPOSE:
! domain partitioning algorithm
      subroutine set_partition_dist
      icdist = 1
      iwdist = 1
      ! associate the local pointers to 'stIcsv'
      call datpol_init_scalar_loc

      ! partition only on grid cells
      !call set_dist_clls
      ! partition on all the blocks
      call set_dist_blocks
      end subroutine set_partition_dist
! PURPOSE:
! distribute the grid cells into each processor
      subroutine set_dist_clls
      use stmapimts, only : metis_api
      use stmgeomet, only : set_loc_clls_size, set_locsize_serial, &
      set_locsize_serial
      if (nprocs == 1) then
          ibdist = 1
          ibcatl = 1
          call set_locsize_serial(stIcsv(7:9), stIcsv(101:103), stIcsv(104:106), &
          stIcsv(107:109), stIcsv(110:112), stIcsv(113:115))
      else
        ! give the graph of grid cells and weighted by transmissibility between grid cells
        call metis_api(0,rank,Ngcll, Ngedges, ivadj, iadjncy, detrans, nprocs, icdist)
        ! 1, metis_recursive_bisection 
        ! 0, metis_kway, keep contiguous of the partition 

      ! based on the partition results, calculate the grid cells type
      ! (inner|border|external) 
      ! count grid cells no. of different type.
        call set_loc_clls_size(deck, Ngcll, icdist, Ngedges, ivadj, iadjncy, detrans, &
        iccatl, Nlcinn, Nlcbrd, Nlcext, Nlcint, Nlncs)
      endif
      end subroutine set_dist_clls
! PURPOSE:
! distribute all the blocks (grid cell + wellblock) into each processor,
! the edge-cut weighted by transmissibility & well index.
!
! CSR-graph of grid (Ngcll, Ngedges, ivadj, iadjncy, detrans)
! CSR-graph of well (Nwidx, Nlyid, ilyadj, icidly, dlywitrs)
      subroutine set_dist_blocks
      use stmapimts, only : metis_api
      use stmgeomet, only : combine_cwll_adj, set_locvert_cate, update_locvert_size, &
      set_locsize_serial
      use stmvtkxml, only : prt_csr
      integer           :: negds
      integer, dimension(Ngidx+1)         :: xadj
      integer, dimension(Ngedges+2*Nlyid)       :: adjncy
      real(STDD), dimension(Ngedges+2*Nlyid)    :: adjwgt

      negds = Ngedges + 2 * Nlyid

      call combine_cwll_adj(Ngcll, Ngedges, ivadj, iadjncy, detrans, &
      Nwidx, Nlyid, ilyadj, icidly, dlywitrs,  &
      Ngidx, negds, xadj, adjncy, adjwgt)
      
      call prt_csr(FUNIT_OUT, Ngidx, negds, xadj, adjncy, adjwgt,"trans-wellindex", 1)
      
      if(nprocs == 1) then
          ! define the local as the global values
          !call set_locsize_serial(stIcsv(7:9), stIcsv(101:103), stIcsv(104:106), &
          !stIcsv(107:109), stIcsv(110:112), stIcsv(113:115))
          ibdist = 1
          ibcatl = 1
      else
        call metis_api(1,rank,Ngidx, negds, xadj, adjncy, adjwgt, nprocs, ibdist)
      ! 1, metis_recursive_bisection 
      ! 0, metis_kway, failed !!!

      ! based on the partition results, calculate the blocks type
      ! (inner|border|external) and count blocks no.
      ! calculate the blocks type (inner|border|external) 
        call set_locvert_cate(deck, Ngidx, ibdist, negds, xadj, adjncy, adjwgt, ibcatl)
      endif

      ! local internal(inner+border) cells no. 
      ! local (internal +external) cells no.
      ! count cells & wells no.
      call update_locvert_size(Ngcll, iccatl, Nlcinn, Nlcbrd, Nlcext, Nlcint, Nlncs) 
      call update_locvert_size(Nwidx, iwcatl, Nlwinn, Nlwbrd, Nlwext, Nlwint, Nlnwls)
      ! count blocks no.
      call update_locvert_size(Ngidx, ibcatl, Nlbinn, Nlbbrd, Nlbext, Nlbint, Nlnbls)

      end subroutine set_dist_blocks
! after get the size of local cells/wells/blocks on current subdeck
      subroutine datpol_init_vect_loc
      szIloc = Nlncs + Nlnwls
      call palloc_i(dpIloc, szIloc, stErrs(15))
      call assoptr_i(dpIloc, szIloc, ofIloc, ibgid, Nlnbls)
      call setptr_i(dpIloc, szIloc, ofIloc, icgid, Nlncs)
      call setptr_i(dpIloc, szIloc, ofIloc, iwgid, Nlnwls)

      end subroutine datpol_init_vect_loc
! after get subdomains number, initialize the vectors relate to communication
      subroutine pmpi_init_vect
      !if master processor
      if(rank /= MASTER) then
          pmnics = 4*nprocs
          call palloc_i(pmIcsv, pmnIcs, stErrs(13))
          call setptr_i(pmIcsv, pmnIcs, pmofIc, pmInec, nprocs)
          call setptr_i(pmIcsv, pmnIcs, pmofIc, pmIdisec, nprocs)
          call setptr_i(pmIcsv, pmnIcs, pmofIc, pmInew, nprocs)
          call setptr_i(pmIcsv, pmnIcs, pmofIc, pmIdisew, nprocs)
        
          ! must associate 'pmInmas' with a memory, otherwise, the following
          ! subroutine pmpi_master_collect_i, will has an dummy parameter
          ! undefine. 
          ! NOTE: here associate it with an arbitrary memory. Assuming it will
          ! not be read and written on NON-MASTER processor
          call assoptr_i(stIgvs, szIdp, 0, pmInmas, 6*nprocs)
          return
      endif
      pmnIcs = nprocs * 10
      call palloc_i( pmIcsv, pmnIcs, stErrs(13))
      
      call setptr_i(pmIcsv, pmnIcs, pmofIc, pmInec, nprocs)
      call setptr_i(pmIcsv, pmnIcs, pmofIc, pmIdisec, nprocs)
      call setptr_i(pmIcsv, pmnIcs, pmofIc, pmInew, nprocs)
      call setptr_i(pmIcsv, pmnIcs, pmofIc, pmIdisew, nprocs)

      call assoptr_i( pmIcsv,  pmnIcs, pmofIc, pmInmas, 6*nprocs)
      call setptr_i( pmIcsv, pmnIcs, pmofIc, pmInupcl, nprocs)
      call setptr_i( pmIcsv, pmnIcs, pmofIc, pmInupwl, nprocs)
      call setptr_i( pmIcsv, pmnIcs, pmofIc, pmInupbk, nprocs)

      call setptr_i( pmIcsv, pmnIcs, pmofIc, pmInlcll, nprocs)
      call setptr_i( pmIcsv, pmnIcs, pmofIc, pmInlwll, nprocs)
      call setptr_i( pmIcsv, pmnIcs, pmofIc, pmInlblk, nprocs)

      end subroutine pmpi_init_vect
! PURPOSE:
! de-associate the pointers
! release the memory
      subroutine datpol_free

      call pmpi_free
      NULLIFY(iblcid, iclcid, iwlcid)

      NULLIFY(ictind, icsecd, icgcat, icdist, iccatl)
      NULLIFY(iwdist, iwdir, iwcatl)
      NULLIFY(ibdist, ibcatl)
      NULLIFY(ixly, iyly, izly, ialy, icidly)
      NULLIFY(ivadj, ilyadj, iadjncy)

      NULLIFY(dcdx, dcdy, dcdz, ds3d)
      NULLIFY(dcareai, dcareaj, dcareak, area3d)
      NULLIFY(dcpor, dcpermi, dcpermj, dcpermk, perm3d)
      NULLIFY(detrans, dethermtrs)
      NULLIFY(dwllrad, dwllgeof, dwllfrac)
      NULLIFY(dlylenfr, dlyskin, dlywitrs, dlywitherm)
      NULLIFY(dchtop, devnts)

      deallocate(stDgvs, STAT = stErrs(11))
      deallocate(stIgvs, STAT = stErrs(12))

      call datpol_free_loc
      end subroutine datpol_free
! PURPOSE:
! mpi related memo.
      subroutine pmpi_free
      if(rank /= MASTER) then
          NULLIFY(pmInec, pmIdisec, pmInew, pmIdisew)
          NULLIFY(pmInmas)
          
          deallocate(dpIloc, STAT=stErrs(16))
          return
      endif
      NULLIFY(pmInec, pmIdisec, pmInew, pmIdisew)

      NULLIFY(pmInmas, pmInlcll, pmInlwll, pmInlblk)
      NULLIFY(pmInupcl, pmInupwl, pmInupbk)
      deallocate(pmIcsv, STAT = stErrs(14))
      end subroutine pmpi_free
      subroutine datpol_free_loc
      NULLIFY(ibgid, icgid, iwgid)
      deallocate(dpIloc, STAT=stErrs(16))
      end subroutine datpol_free_loc
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
