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
      else
          call no_adjclls_five(nx, ny, nz, adjclls)
      endif
      num = sum(adjclls)
      end subroutine no_totl_connections
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
      
      end module stmgeomet
