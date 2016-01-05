      module stmgeomet
      use stmheader
      contains
      
! local variables equals to global variables if no domain partitioning
      elemental subroutine set_locsize_serial(totln, ninn, nbrd, next, ninternal, nloc)
      integer, intent(in)                   :: totln
      integer, intent(out)                  :: ninn, nbrd, next
      integer, intent(out)                  :: ninternal, nloc
      nloc      = totln   ! local grid cells number = total
      ninternal = totln
      ninn      = totln ! inner
      nbrd      = 0 ! border
      next      = 0 ! external
      end subroutine set_locsize_serial
! count the local inner, border, external grid cells no.
      subroutine set_loc_clls_size(deck, nvtxs, ddist, nadys, iadj, adjncy, &
      adjwgt, ivcatl, ninn, nbrd, next,ninternal, nloc)
      integer, intent(in)                   :: deck
      integer, intent(in)                   :: nvtxs
      integer, dimension(nvtxs), intent(in) :: ddist
      integer, intent(in)                   :: nadys
      integer, dimension(nvtxs+1), intent(in)   :: iadj
      integer, dimension(nadys), intent(in)     :: adjncy
      real(STDD), dimension(nadys), intent(in)  :: adjwgt
      integer, dimension(nvtxs), intent(out)    :: ivcatl
      integer, intent(out)                  :: ninn, nbrd, next
      integer, intent(out)                  :: ninternal, nloc
      call set_locvert_cate(deck, nvtxs, ddist, nadys, iadj, adjncy, &
      adjwgt, ivcatl)
      call update_locvert_size(nvtxs, ivcatl, ninn, nbrd, next,      &
      ninternal, nloc)
      end subroutine set_loc_clls_size
! counts no.
      subroutine update_locvert_size(nvtxs, ivcatl, ninn, nbrd, next, & 
      ninternal, nloc)
      integer, intent(in)                   :: nvtxs
      integer, dimension(nvtxs), intent(in) :: ivcatl
      integer, intent(out)                  :: ninn, nbrd, next
      integer, intent(out)                  :: ninternal, nloc
      call set_locvert_size(nvtxs, ivcatl, ninn, nbrd, next)
      ninternal = ninn + nbrd
      nloc      = ninternal + next
      end subroutine update_locvert_size
! set up local index by the sequence, & external vertex is re-sorted by
! its subdeck
! input:
! nvg, global vertex no.
! ivcatl, vertex category based on local subdeck
! dist,   vertex distributed subdeck
! ninn,   local inner vertex no.
! nbrd,   local border vertex no.
! nvl,  local vertex no.
! output:
! ivlid(nvg): the local index of a given global order (natural order)
! ivgid(nvl): the natural index of a given local order (inner/border/external)
      subroutine set_locid(nvg,ivcatl, dist, ninn, nbrd, next,      &
      ivlid, nvl, ivgid)
      integer, intent(in)                   :: nvg
      integer, dimension(nvg), intent(in)   :: ivcatl
      integer, dimension(nvg), intent(in)   :: dist
      integer, intent(in)                   :: ninn, nbrd, next
      integer, dimension(nvg), intent(out)  :: ivlid
      integer, intent(in)                   :: nvl
      integer, dimension(nvl), intent(out)  :: ivgid

      integer :: off(4)
      call nums2displs_i(4, (/ninn, nbrd, next, 0/), off)
      call set_locid_raw(nvg,ivcatl, off(1:3), ivlid, nvl, ivgid)
      if(next>0) call sort_locid_ext(next, ivgid(off(3)+1: off(4)), nvg, dist, ivlid)
      end subroutine set_locid
! set up the local index by the sequence (inner, border, external) 
      subroutine set_locid_raw(nvg, ivcatl, disp, ivlid, nvl, ivgid)
      integer, intent(in)                   :: nvg
      integer, dimension(nvg), intent(in)   :: ivcatl
      integer, dimension(3), intent(in)     :: disp
      integer, dimension(nvg), intent(out)   :: ivlid
      integer, intent(in)                   :: nvl
      integer, dimension(nvl), intent(out)   :: ivgid
      integer :: h, off(3)
      off = disp

      do h = 1, nvg
        select case(ivcatl(h))
        case (1) ! inner
            off(1)          = off(1) + 1
            ivlid(h)        = off(1)
            ivgid(off(1))   = h
        case (2) ! border
            off(2)          = off(2) + 1
            ivlid(h)        = off(2)
            ivgid(off(2))   = h
        case (3) ! external
            off(3)          = off(3) + 1
            ivlid(h)        = off(3)
            ivgid(off(3))   = h
        case (0) ! not on current subdeck
            ivlid(h)        = 0
        end select
      enddo
      end subroutine set_locid_raw
! sort the external vertex sequence according to where(which subdeck) it comes from.
! then, the external blocks are stored contiguously
! insertion sort
! note: ivgid(next) only the external parts are passed in & out
! ivgid(next) stores the natural index of external vertex subdeck by subdeck
      subroutine sort_locid_ext(next, ivgid, nvg, dist, ivlid)
      integer, intent(in)                       :: next
      integer, dimension(next), intent(inout)   :: ivgid
      integer, intent(in)                       :: nvg
      integer, dimension(nvg), intent(in)       :: dist
      integer, dimension(nvg), intent(inout)    :: ivlid

      integer :: i,j,iptr
      integer :: temp
      do i = 1, next - 1
        iptr = i
        do j = i + 1, next
            if( dist(ivgid(j)) < dist(ivgid(iptr)) ) then
                iptr = j

                temp            = ivlid( ivgid(i) )
                ivlid(ivgid(i)) = ivlid(ivgid(iptr))
                ivlid(ivgid(iptr)) = temp

                temp        = ivgid(i)
                ivgid(i)    = ivgid(iptr)
                ivgid(iptr) = temp
                
                exit
            endif ! swap
        enddo ! j
      enddo ! i
      end subroutine sort_locid_ext
! external vertex number and displs from different adjacent subdecks
      subroutine set_numdisp_extofadjdeck(nvg, dist, ivcatl, ndeck, num, displ)
      integer, intent(in)                   :: nvg
      integer, dimension(nvg), intent(in)   :: dist
      integer, dimension(nvg), intent(in)   :: ivcatl
      integer, intent(in)                   :: ndeck
      integer, dimension(ndeck), intent(out):: num
      integer, dimension(ndeck), intent(out):: displ
      call set_num_extofadjdeck(nvg,dist,ivcatl, ndeck, num)
      call nums2displs_i(ndeck,num,displ)
      end subroutine set_numdisp_extofadjdeck
! cout the external vertex number from different adjacent subdecks
      subroutine set_num_extofadjdeck(nvg, dist, ivcatl, ndeck, numext)
      integer, intent(in)                   :: nvg
      integer, dimension(nvg), intent(in)   :: dist
      integer, dimension(nvg), intent(in)   :: ivcatl
      integer, intent(in)                   :: ndeck
      integer, dimension(ndeck), intent(out):: numext

      integer :: h
      numext = 0
      forall(h=1:nvg, ivcatl(h) == 3) 
        numext(dist(h)) = numext(dist(h)) + 1
      end forall
      end subroutine set_num_extofadjdeck
! count the local inner, border, external grid cells no.
      subroutine set_locvert_size(nvtxs, ivcatl, ninn, nbrd, next)
      integer, intent(in)                   :: nvtxs
      integer, dimension(nvtxs), intent(in) :: ivcatl
      integer, intent(out)                  :: ninn, nbrd, next

      integer  :: h, iarr(3)
      iarr = 0
!$OMP PARALLEL if(nvtxs>ST_OMP_LIMIT)
!$OMP DO REDUCTION(+:iarr)
      do h = 1, nvtxs
        select case( ivcatl(h))
        case (1) ! inner
            iarr(1) = iarr(1) + 1
        case (2) ! border
            iarr(2) = iarr(2) + 1
        case (3) ! external
            iarr(3) = iarr(3) + 1
        end select
      enddo
!$OMP END DO
!$OMP END PARALLEL
      ninn = iarr(1) 
      nbrd = iarr(2) 
      next = iarr(3)
      end subroutine set_locvert_size
! calculate the block category based on local subdeck
!   * 0, not on current subdeck
!   * 1, inner 
!   * 2, border
!   * 3, external
!   * 1:2 internal
!   * 1:3 local
      subroutine set_locvert_cate(deck, nvtxs, ddist, nadys,        &
      iadj, adjncy, adjwgt, ivcatl)
      integer, intent(in)                   :: deck
      integer, intent(in)                   :: nvtxs
      integer, dimension(nvtxs), intent(in) :: ddist
      integer, intent(in)                   :: nadys
      integer, dimension(nvtxs+1), intent(in)   :: iadj
      integer, dimension(nadys), intent(in)     :: adjncy
      real(STDD), dimension(nadys), intent(in)  :: adjwgt
      integer, dimension(nvtxs), intent(out)    :: ivcatl

      integer :: h, n
      logical :: border
      ivcatl = 0
!$OMP PARALLEL PRIVATE(n, border) if(nvtxs>ST_OMP_LIMIT)
!$OMP FLUSH(ddist, iadj, adjncy, adjwgt)
!$OMP DO
      do h = 1, nvtxs
        border = .false.
        if( ddist(h) /= deck ) cycle
        do n = iadj(h), iadj(h+1) - 1
          if( ddist(adjncy(n)) /= deck) then ! external
            if( adjwgt(n) > ST_TOL_EPSILON ) ivcatl( adjncy(n) ) = 3  ! tranmissibility > 0
            border = .true.
          endif
        enddo
        if(border) then
            ivcatl(h) = 2
        else
            ivcatl(h) = 1
        endif
      enddo
!$OMP END DO
!$OMP FLUSH(ivcatl)
!$OMP END PARALLEL
      end subroutine set_locvert_cate
! calculate the global connections number in the whole region
! (1:nx,1:ny,1:nz)
      subroutine no_totl_connections(nx,ny,nz, forn, num)
      integer, intent(in)   :: nx, ny, nz
      integer, intent(in)   :: forn ! five or nine stencil
      integer, intent(out)  :: num ! total connection numbers

      integer, dimension(nz*ny*nx) :: adjclls
      if( forn == 1) then
          ! todo
      else
          call no_adjclls_five(nx, ny, nz, adjclls)
      endif
      num = sum(adjclls)
      end subroutine no_totl_connections
! calculate the total connections number 'nadys' and index array 'iadj' into CSR
! adjncy(nadys) 
      subroutine set_totl_connect_size(nx,ny,nz,forn, nadys, iadj)
      integer, intent(in)   :: nx, ny, nz
      integer, intent(in)   :: forn ! five or nine stencil
      integer, intent(out)  :: nadys
      integer, dimension(nz*ny*nx+1), intent(out)  :: iadj

      integer, dimension(nz*ny*nx) :: adjclls
      integer :: n
      n = nz * ny * nx
      if( forn == 1) then
          ! todo
      else
          call no_adjclls_five(nx, ny, nz, adjclls)
      endif
      call nadjcny2iadj(n,adjclls,nadys,iadj)
      end subroutine set_totl_connect_size
! calculate the CSR-graph 'adjncy(nadys)' and corresponding
! tranmissibility 'adjtrans(nadys)'
      subroutine set_adjncy_trans(nx,ny,nz,forn, nadys, iadj, perm, area, step, &
      adjncy, adjtrans)
      integer, intent(in)   :: nx, ny, nz
      integer, intent(in)   :: forn ! five or nine stencil
      integer, intent(in)   :: nadys
      integer, dimension(nz*ny*nx+1), intent(in)  :: iadj
      real(STDD), dimension(3*nz*ny*nx), intent(in) :: perm, area, step
      integer, dimension(nadys), intent(out)      :: adjncy
      real(STDD), dimension(nadys), intent(out)      :: adjtrans

      integer :: nclls
      integer, dimension(nadys)         :: adjdir
      nclls = nz * ny * nx
      if( forn == 1) then
          ! todo
      else
          call adjncy_clls_five(nx,ny,nz, nadys, iadj, adjncy, adjdir)
          call adjncy_clls_trans_five_cart(nclls, nadys, iadj, adjncy, &
          adjdir, perm, area, step, adjtrans)
      endif
      end subroutine set_adjncy_trans
! calculate the index into adjcny(nadys) that is the begining of the
! adjacency list of vertices
! input:
! adjclls(nvtxs) : adjacency blocks number of each block
! output:
! nadys, total edges|connections no. ( no. is count repeatedly since adjclls(nvtxs)
! stores for both points )
! iadj(nvtxs+1): the beginning index of csr adjacncy list
      subroutine nadjcny2iadj(nvtxs, adjclls, nadys, iadj)
      integer, intent(in)                       :: nvtxs
      integer, dimension(nvtxs), intent(in)     :: adjclls
      integer, intent(out)                      :: nadys
      integer, dimension(nvtxs+1), intent(out)  :: iadj

      integer :: i
      iadj(1) = 1
      do i = 1, nvtxs
        iadj(i+1) = iadj(i) + adjclls(i)
      enddo
      nadys = iadj(nvtxs+1) - 1  !! == sum(adjclls)
      end subroutine nadjcny2iadj
! count the adjacent grid-cells no. for each cell at five-stencial
! input:
! nx, the grid cells no. along I/x direction
! ny,                          J/y
! nz,                          K/z
! The natural order (label starts from x, then y, following z)
! output:
! nadjc(nz*ny*nx) : the adjacent cells no. for each cell (in natural order) 
      subroutine no_adjclls_five(nx, ny, nz, nadjc)
      integer, intent(in)                       :: nx, ny, nz
      integer, dimension(nz*ny*nx), intent(out) :: nadjc

      integer       :: i, j, k, zf, yf, id
      do k = 1, nz
        zf = (k-1) * ny * nx
        do j = 1, ny
            yf = (j-1) * nx
            do i = 1, nx
                id = zf + yf + i 
                nadjc(id) = no_adjclls_five_ofcell( i, 1, nx, j, 1, ny, k, 1, nz)
            enddo
        enddo
      enddo
      end subroutine no_adjclls_five
! calculate the adjacency lists of the grid cells for five-stencial FD
! the adjacency direction lists are stored in 'adjdir', with acceptable
! value: 
!   0, I direction
!   1, J direction
!   2, K direction
      subroutine adjncy_clls_five(nx,ny,nz, nadys, iadj, adjncy, adjdir)
      integer, intent(in)                               :: nx, ny, nz
      integer, intent(in)                               :: nadys
      integer, dimension(nz*ny*nx+1), intent(in)        :: iadj
      integer, dimension(nadys), intent(out)            :: adjncy
      integer, dimension(nadys), intent(out)            :: adjdir

      integer           :: i,j,k, gid, n, a0, a1
      do k = 1, nz
        do j = 1, ny
            do i = 1, nx
            gid = (k-1) * nx * ny + (j-1) * nx + i
            a0 = iadj(gid)
            a1 = iadj(gid+1)
            n = a1 - a0
            call adjncy_five_ofcell(i, 1, nx, j, 1, ny, k, 1, nz, n,  &
            adjncy(a0 : a1-1), adjdir(a0 : a1-1) )
            enddo
        enddo
      enddo
      end subroutine adjncy_clls_five
! calculate the adjacency cells and directions lists of a given grid
! cell with UID(ui,uj,uk) /in ( (i0,ni), (j0,nj), (k0, nk) )
! input:
! nadys, adjacency cells number of given grid cell
! output:
! adjncy(nadys)  : adjacency cells gid in natural order
! adjdir(nadys)  : direction of <h,n>
      pure subroutine adjncy_five_ofcell(ui, i0, ni, uj, j0, nj, uk, k0, &
      nk, nadys, adjncy, adjdir)
      integer, intent(in)           :: ui, uj, uk
      integer, intent(in)           :: i0, j0, k0
      integer, intent(in)           :: ni, nj, nk
      integer, intent(in)           ::  nadys
      integer, dimension(nadys), intent(out)        :: adjncy
      integer, dimension(nadys), intent(out)        :: adjdir
      integer :: off
      off = 1
      if ( uk > k0 ) call adjncy_five_ofcell_ofdim(ui,uj,uk-1,ni,nj,2,adjncy(off),adjdir(off),off)
      if ( uj > j0 ) call adjncy_five_ofcell_ofdim(ui,uj-1,uk,ni,nj,1,adjncy(off),adjdir(off),off)
      if ( ui > i0 ) call adjncy_five_ofcell_ofdim(ui-1,uj,uk,ni,nj,0,adjncy(off),adjdir(off),off)
      if ( ui < ni ) call adjncy_five_ofcell_ofdim(ui+1,uj,uk,ni,nj,0,adjncy(off),adjdir(off),off)
      if ( uj < nj ) call adjncy_five_ofcell_ofdim(ui,uj+1,uk,ni,nj,1,adjncy(off),adjdir(off),off)
      if ( uk < nk ) call adjncy_five_ofcell_ofdim(ui,uj,uk+1,ni,nj,2,adjncy(off),adjdir(off),off)

      end subroutine adjncy_five_ofcell
      pure subroutine adjncy_five_ofcell_ofdim( ui, uj, uk, nx, ny, idir, adjncy, adjdir, off)
      integer, intent(in)               :: ui,uj, uk
      integer, intent(in)               :: nx,ny
      integer, intent(in)               :: idir
      integer, intent(out)              :: adjncy
      integer, intent(out)              :: adjdir
      integer, intent(inout)            :: off
      adjdir = idir
      adjncy = U2GID(ui,uj,uk,nx,ny)
      off = off + 1
      end subroutine adjncy_five_ofcell_ofdim
! calculate the geometry category of the grid cell in Bits info
! INPUT:
! nx,ny,nz: TOTAL GLOBAL grid cells number along x, y, z direction
! output:
! BTEST(gcat, 0) 
! BTEST(gcat, 1) I- geometric boundary
! BTEST(gcat, 2) I+
! BTEST(gcat, 3) J-
! BTEST(gcat, 4) J+
! BTEST(gcat, 5) under-burden, (kdir up: K-, down: K+)
! BTEST(gcat, 6) over-burden (kdir up: K+, down: K-)
      subroutine set_clls_gcate(nx,ny,nz, gcat)
      integer, intent(in)                           :: nx, ny, nz
      integer, dimension(nz*ny*nx), intent(out)     :: gcat

      integer :: n, i, s, e, k
      n = nz * ny * nx
      do i = 1, n, nx
        gcat(i) = IBSET(gcat(i), 1)
        gcat(i+nx-1) = IBSET(gcat(i+nx-1), 2)
      enddo
      ! K = 0 layer
      s = ny * nx
      gcat(1 : s) = IBSET(gcat(1 : s), 6) ! kdir down
      ! K = nz layer
      e = (nz-1)*s + 1 
      gcat(e : n) = IBSET(gcat(e: n), 5)
      ! switch if kdir up

      ! J
      do k = 1, nz
        e = (k-1) * s
        gcat( e+1 : e+nx) = IBSET(gcat(e+1 : e+nx), 3)
        e = k * s - nx
        gcat( e+1 : e+nx) = IBSET(gcat(e+1 : e+nx), 4)
      enddo

      end subroutine set_clls_gcate
! calculate the transmissibility AK/h of each connection. 
! input: 
!      The graph info. CSR (nclls, nadys, iadj, adjncy) format 
!       properties:
!       perm = (permi, permj, permk) 
!       area = (areai, areaj, areak) 
!   the area perpendicular to the direction
!       step = (dx, dy, dz)
! output:
!   transmissibility AK/h
      subroutine adjncy_clls_trans_five_cart(nclls, nadys, iadj, adjncy, adjdir, &
      perm, area, step, adjtrans)
      integer, intent(in)                               :: nclls
      integer, intent(in)                               :: nadys
      integer, dimension(nclls+1), intent(in)           :: iadj
      integer, dimension(nadys), intent(in)             :: adjncy
      integer, dimension(nadys), intent(in)             :: adjdir
      real(STDD), dimension(3*nclls), intent(in)        :: perm
      real(STDD), dimension(3*nclls), intent(in)        :: area
      real(STDD), dimension(3*nclls), intent(in)        :: step
      real(STDD), dimension(nadys), intent(out)         :: adjtrans

      integer       :: i, l, j, dir
      real(STDD)    :: permh, permn, areah, arean, steph, stepn

      do i = 1, nclls
        do l = iadj(i), iadj(i+1) - 1
            j = adjncy(l)
            dir = adjdir(l)
            ! <i,j>
            call set_para_trans_cart_five(dir, i, j, nclls, perm, area, step, &
            permh, permn, areah, arean, steph, stepn)
            call harm_cart_trans_five_hete(permh, permn, areah, arean, &
            steph, stepn, adjtrans(l))
        enddo
      enddo 
      end subroutine adjncy_clls_trans_five_cart
! assign the corresponding permeability, area, and distance of an edge
! of grid. 
! input:
!   <h,n>, index of the two vertices of the edge in natural order
!  dir, which direction of the edges along with. 
!   0, I direction
!   1, J direction
!   2, K direction
! perm[], area[], step[], store the permeability, area, and distance of
! the cartisian grids, with continous memory in I,J,K directions. i.e.
! perm = /permi, permj, permk/
! output:
! corresponding perm, area, step to compute the transmissibility of
! connection <h,n>
      subroutine set_para_trans_cart_five(dir, h, n, nclls, perm, area, &
      step, permh, permn, areah, arean, steph, stepn)
      integer, intent(in)       :: dir
      integer, intent(in)       :: h, n ! natural order of home|neigh clls
      integer, intent(in)       :: nclls
      real(STDD), dimension(3*nclls) , intent(in)       :: perm
      real(STDD), dimension(3*nclls) , intent(in)       :: area
      real(STDD), dimension(3*nclls) , intent(in)       :: step
      real(STDD), intent(out)                           :: permh, permn
      real(STDD), intent(out)                           :: areah, arean
      real(STDD), intent(out)                           :: steph, stepn
      permh = perm( dir * nclls + h)
      areah = area( dir * nclls + h)
      steph = step( dir * nclls + h)
      permn = perm( dir * nclls + n)
      arean = area( dir * nclls + n)
      stepn = step( dir * nclls + n)
      end subroutine set_para_trans_cart_five

! count neighbors of the given grid-cells@ UID(ui,uj,uk), global
! region((i0,ni),(j0, nj), (k0, nk)
      pure integer function no_adjclls_five_ofcell(ui, i0, ni,uj, j0, nj, uk, k0, nk) result(num)
      integer, intent(in)              :: ui, uj, uk
      integer, intent(in)              :: i0, j0, k0
      integer, intent(in)              :: ni, nj, nk
      num = 0
      num = num + no_adjclls_five_ofcell_ofdim(ui, i0, ni)
      num = num + no_adjclls_five_ofcell_ofdim(uj, j0, nj)
      num = num + no_adjclls_five_ofcell_ofdim(uk, k0, nk)
      end function no_adjclls_five_ofcell 
! count the neighbors no. of fixed index(ui) at given dimension (i0<= ui <=ni)
      pure integer function no_adjclls_five_ofcell_ofdim(ui, i0, ni) result(neigh)
      integer, intent(in)              :: ui, i0, ni
      neigh = 0
      if( ui < ni ) then ! ui<ni, normal case
        if( ui > i0) then ! inner cells
          neigh = 2
        else if(ui == i0) then  !  lower bouindary
          neigh = 1
        else
          return
        endif
      else if(ui == ni) then ! ui>ni
        if(i0 == ni) then ! this dimension does not exist
          return
        else
          neigh = 1
        endif
      else ! ui > ni
        return
      endif
      return
      end function no_adjclls_five_ofcell_ofdim
! calculate the transmissibility between two grid-cells
! input: 
! perm_h, perm_n: permeabilities in current(h) and neighboring(n) cells
! area_h, area_n: area perpendicular to the flow direction
! step_h, step_n: grid stride (step size) along the flow direction
! output:
!   transmissibility, units = unit[perm] * unit[step]
      subroutine harm_cart_trans_five_hete(permh, permn, areah, arean, steph, stepn, trans)
      real(STDD), intent(in)    :: permh, permn 
      real(STDD), intent(in)    :: areah, arean
      real(STDD), intent(in)    :: steph, stepn
      real(STDD), intent(out)   :: trans

      real(STDD)                :: fh, fn
      fh = areah * permh 
      fn = arean * permn 
      trans = 2 * fh * fn / (fh * stepn + fn * steph)
      ! CARTISIAN GRID
      ! The transmissibility at the interface <h, n> (h & n are cell index)
      ! is the average 
      ! \frac{1}{\trans_{<h,n>}}  = \left(\frac{1}{\trans_h}
      ! + \frac{1}{\trans_n} \right)
      ! \trans = \frac{area * perm}{0.5 * step}
      end subroutine harm_cart_trans_five_hete
! calculate the area = length x width
      elemental subroutine set_area(length, width, area)
      real(STDD), intent(in)            :: length, width
      real(STDD), intent(out)           :: area
      area = length * width
      end subroutine set_area
! calculate the center of the grid cell (Cartesian)
! ni = NX, nj = NY, nk = NZ
      subroutine set_center(ni, nj, nk, step, d, center)
      integer, intent(in)                               :: ni,nj,nk
      real(STDD), dimension(ni,nj,nk), intent(in)       :: step
      integer, intent(in)                               :: d
      real(STDD), dimension(ni,nj,nk), intent(out)      :: center

      real(STDD), dimension(ni*nj*nk)                   :: temp
      real(STDD), dimension(:), allocatable             :: dstep, dcet
      integer :: ndim(3), adim(3)
      ndim       = (/ni,nj,nk/)
      adim       = CSHIFT( ndim, (d-1) )
      ! d=1, adim = (NX,NY,NZ)
      !  =2, adim = (NY,NZ,NX)
      !  =3, adim = (NZ,NX,NY)
      ! ref(nr) = ref(ndim(2), ndim(3)) 
      allocate( dstep(ndim(d)), dcet(ndim(d)) )

      select case (d)
      case (1)
          dstep = step( :, 1, 1)
      case (2)
          dstep = step( 1, :, 1)
      case (3)
          dstep = step( 1, 1, :)
      end select

      call step2center_d(ndim(d), dstep, dcet)      
      call vector2array_d(d, ndim(d), dcet, adim(2), adim(3), temp)
      center = RESHAPE(temp, SHAPE(center))
      deallocate(dstep, dcet)
      end subroutine set_center
      subroutine set_center_htop(ni, nj, nk, step, ref, center)
      integer, intent(in)                               :: ni,nj,nk
      real(STDD), dimension(ni,nj,nk), intent(in)       :: step
      real(STDD), dimension(ni,nj), intent(in)          :: ref
      real(STDD), dimension(ni,nj,nk), intent(out)      :: center
      integer :: i,j
      call set_center(ni,nj,nk,step, 3, center)
      forall(i =1:ni, j=1:nj)
        center(i,j,:) = center(i,j,:) + ref(i,j)
      end forall
      end subroutine set_center_htop
! calculate the cell center with cells dimensional step 
      subroutine step2center_d(nd, step, center)
      integer, intent(in)                               :: nd
      real(STDD), dimension(nd), intent(in)             :: step
      real(STDD), dimension(nd), intent(out)            :: center
      real(STDD), dimension(nd) :: disp
      call nums2displs_d(nd, step, disp)
      center = disp + 0.5 * step
      end subroutine step2center_d
! expands vector to be array, ivar|jvar|kvar
! input:
! if d=1, ni = NX, nj=NY, nk=NZ, 
! if d=2, ni = NY, nj=NZ, nk=NX, 
! if d=3, ni = NZ, nj=NX, nk=NY, 
      subroutine vector2array_d(d, ni, vect, nj, nk, arr)
      integer, intent(in)                               :: d
      integer, intent(in)                               :: ni, nj, nk
      real(STDD), dimension(ni), intent(in)             :: vect
      real(STDD), dimension(ni*nj*nk), intent(out)      :: arr

      integer :: ndim(3), d1, d2
      real(STDD), dimension(:,:,:), allocatable     :: temp
      ndim = CSHIFT( (/ni,nj,nk/), (4-d) )
      allocate( temp(ndim(1), ndim(2), ndim(3)) )
      select case (d)
        case (1) ! ivar
            forall(d1 =1: ndim(2), d2=1: ndim(3) )
                temp(:,d1,d2) = vect(:)
            end forall
        case (2) ! jvar
            forall(d1 = 1:ndim(1), d2=1:ndim(3))
                temp(d1, :, d2) = vect(:)
            end forall
        case (3) ! kvar
            forall (d1=1: ndim(1), d2=1:ndim(2))
                temp(d1, d2, :) = vect(:)
            end forall
      end select
      arr = RESHAPE(temp, (/ni*nj*nk/) )
      deallocate(temp)
      end subroutine vector2array_d
! PURPOSE:
! combine the connection lists between grid cells, and between wellbore and
! perforated layers into a total CSR-graph 
! input:
!   CSR-graph of grid cells:
!     ncls, nnzcl, icadj, icncy, dcwgt
!   CSR-graph of well perforated info.
!     nwls, nnzwl, iwadj, iwncy, dwwgt
      subroutine combine_cwll_adj(ncls, nnzcl, icadj, icncy, dcwgt, &
      nwls, nnzwl, iwadj, iwncy, dwwgt, n, nnz, xadj, adjncy, adjwgt)
      ! CSR-graph on grid cells
      integer, intent(in)                       :: ncls, nnzcl
      integer, dimension(ncls+1), intent(in)    :: icadj
      integer, dimension(nnzcl), intent(in)     :: icncy
      real(STDD), dimension(nnzcl), intent(in)  :: dcwgt
      ! CSR-graph on well-grid
      integer, intent(in)                       :: nwls, nnzwl
      integer, dimension(nwls+1), intent(in)    :: iwadj
      integer, dimension(nnzwl), intent(in)     :: iwncy
      real(STDD), dimension(nnzwl), intent(in)  :: dwwgt

      integer, intent(in)                       :: n, nnz
      integer, dimension(n+1) , intent(out)     :: xadj
      integer, dimension(nnz) , intent(out)     :: adjncy
      real(STDD), dimension(nnz) , intent(out)  :: adjwgt

      integer, dimension(ncls+1)                :: icwadj
      integer, dimension(nnzwl)                 :: icwncy
      real(STDD), dimension(nnzwl)              :: dcwwgt

      integer, dimension(ncls) ::n1adj, n2adj, nadj
      integer :: i
      
      call csr2csc(nr=nwls, nc=ncls, nnz=nnzwl, ap=iwadj, aj=iwncy,  &
      av=dwwgt, bp=icwadj, bi=icwncy, bv=dcwwgt)

      n1adj = icadj(2 :ncls+1) - icadj(1 :ncls)
      n2adj = icwadj(2:ncls+1) - icwadj(1:ncls)
      nadj = n1adj + n2adj

      xadj = 0
      xadj(1) = 1
      do i = 1, ncls 
        xadj(i+1) = xadj(i) + nadj(i)
        
        adjncy( xadj(i) : xadj(i) + n1adj(i)-1 ) = icncy(icadj(i)   : icadj(i+1)-1)
        ! order the wellbore behind all grid cells
        adjncy( xadj(i)+n1adj(i) : xadj(i+1)-1 ) = icwncy(icwadj(i) : icwadj(i+1)-1) + ncls

        adjwgt( xadj(i) : xadj(i) + n1adj(i)-1 ) = dcwgt(icadj(i)   : icadj(i+1)-1)
        adjwgt( xadj(i)+n1adj(i) : xadj(i+1)-1 ) = dcwwgt(icwadj(i) : icwadj(i+1) - 1)
      enddo
      ! order the wellbore behind all grid cells
      do i = 1, nwls
        xadj(ncls+1 + i) = xadj(ncls + i) + iwadj(i+1) - iwadj(i)

        adjncy(xadj(ncls+i) : xadj(ncls+1+i)-1) = iwncy(iwadj(i) : iwadj(i+1)-1)

        adjwgt(xadj(ncls+i) : xadj(ncls+1+i)-1) = dwwgt(iwadj(i) : iwadj(i+1)-1)
      enddo

      end subroutine combine_cwll_adj
! PURPOSE:
! convert between gid (natural order x-y-z) and UID
      pure integer function U2GID(ui, uj, uk, nx, ny)
      integer, intent(in)       :: ui, uj, uk
      integer, intent(in)       :: nx, ny
      call U2GID_(ui,uj,uk,nx,ny, U2GID)
      end function U2GID
      pure subroutine G2UID(gid, nx, ny, ui, uj, uk)
      integer, intent(in)       :: gid
      integer, intent(in)       :: nx, ny
      integer, intent(out)      :: ui,uj,uk
      integer :: v1, v2
      v1 = gid - 1
      v2 = MOD(v1, nx*ny)
      ui = MOD(v2, nx) + 1
      uj = v2 / nx + 1
      uk = v1 / (nx*ny) + 1
      end subroutine G2UID
      pure subroutine U2GID_(ui,uj,uk,nx, ny, gid)
      integer, intent(in)       :: ui, uj, uk
      integer, intent(in)       :: nx, ny
      integer, intent(out)      :: gid
      gid = (uk-1) * ny * nx + (uj-1) * nx + ui
      end subroutine U2GID_ 
! PURPOSE:
! circulate the direction
! input: 
! i 
! output: 
! j, k
! <i,j,k> = (0, 1, 2),
!           (1, 2, 0)
!           (2, 0, 1)
      subroutine get_perpend_dirs(i,j,k)
      integer, intent(in)           :: i
      integer, intent(out)          :: j,k
      select case (i)
      case (0)
          j = 1
          k = 2
      case (1)
          j = 2
          k = 0
      case (2)
          j = 0
          k = 1
      case default
          j = -1
          k = -1
      end select
      end subroutine get_perpend_dirs
! PURPOSE:
! extract parts of the properties from an array
! mid(m) : positions in arr
      subroutine get_geo_prop_d(n, arr, m, mid, extract)
      integer, intent(in)                   :: n
      real(STDD), dimension(n), intent(in)  :: arr
      integer, intent(in)                   :: m
      integer, dimension(m), intent(in)     :: mid
      real(STDD), dimension(m), intent(out) :: extract
      extract = arr(mid) 
      end subroutine get_geo_prop_d

      end module stmgeomet
