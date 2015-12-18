      module stmgeomet
      use stmheader
      contains
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
! nadys, total edges|connections no. count repeated since adjclls(nvtxs)
! stores for both points
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
      nadys = iadj(nvtxs+1)
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
! BTEST(gcat, 5) K- under-burden
! BTEST(gcat, 6) K+ over-burden
      subroutine set_clls_gcate(nx,ny,nz, gcat)
      integer, intent(in)                           :: nx, ny, nz
      integer, dimension(nz*ny*nx), intent(out)     :: gcat

      integer :: n, i
      n = nz * ny * nx
      do i = 1, n, nx
        gcat(i) = IBSET(gcat(i), 1)
        gcat(i+nx-1) = IBSET(gcat(i+nx-1), 2)
      enddo
      ! todo
      ! J
      ! K
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
      end module stmgeomet
