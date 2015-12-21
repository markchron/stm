      module stmwllbor
      use stmheader
      real(STDD), parameter :: ST_FUCONV = 1.127e-3 
      ! field units, permeability in 'mD', viscosity in 'cp', pressure
      ! in 'psi', length in 'ft', area in 'sq-ft', volumetric rate in
      ! 'stb/d'
      ! mass in 'lb', time in 'hr' density in 'lb/cuft', velocity in
      ! 'ft/sec'
      ! Darcy: perm = 1 Darcy, when u = 1cm/s, \mu = 1 cp and dp/dl = 1 atm/cm
      ! Darcy units, permeability in 'D', viscosity in 'cp', pressure in
      ! 'atm', length in 'cm', velocity 'cm/s'
      ! area in 'cm^2', volumetric rate in 'cc/s', mass in 'gm', time in 'sec', 
      ! density 'gm/cc'
      ! SI units, permeability 'm^2', viscosity in 'kg/m-s', pressure
      ! 'pa', length 'm', velocity 'm/s'
      ! area in 'm^2', volumetric rate 'm3/s', mass in 'kg', time in
      ! 's', density in 'kg/m3'
      ! THE CGS and SI units are impracticably large for the majority of
      ! reservoir rock, the field or Darcy units was devised in which
      ! the permeability would have a more convenient numerical size.
      contains
! calculate the UID (in natural order) of perforated layer of all wells
      subroutine set_wll_cid(nlays, xlay, ylay, zlay, ulay, nx, ny)
      use stmgeomet, only : U2GID
      integer, intent(in)           :: nlays
      integer, dimension(nlays), intent(in) :: xlay, ylay, zlay
      integer, intent(in)           :: nx, ny
      integer, dimension(nlays), intent(out) :: ulay
      
      integer :: i
      forall(i = 1: nlays) 
        ulay(i) = U2GID(xlay(i), ylay(i), zlay(i), nx, ny)
      end forall
      end subroutine set_wll_cid
! calculate equivalent well radius Five-stencial
      elemental subroutine vertical_wll_re_hete( geofac, permi, permj, hi, hj, re)
      real(STDD), intent(in)    :: permi, permj
      real(STDD), intent(in)    :: hi, hj
      real(STDD), intent(in)    :: geofac
      real(STDD), intent(out)   :: re
      real(STDD), parameter     :: para = 2.d0 / sqrt(ST_PI)
      real(STDD)                :: fi, fj
      fi = sqrt(permi/permj) 
      fj = sqrt(permj/permi)

      re = para*geofac * sqrt(fi*hj*hj + fj*hi*hi) / (sqrt(fi)+sqrt(fj))
      end subroutine vertical_wll_re_hete
! calculate the permeability perpendicular to the well direction
      elemental subroutine vertical_wll_perme_hete(permi, permj, perm)
      real(STDD), intent(in)    :: permi, permj
      real(STDD), intent(out)   :: perm
      perm = sqrt(permi * permj)
      end subroutine vertical_wll_perme_hete
! calculate the well-index
! input:
!   wfrac, well fraction (fraction)
!   perm, estimate formation permeability perpendicular to the well direction 
!   length, well length 
!   lthfra, perforated length fraction
!   re, equivalent well radius
!   rw, wellbore radius
!   skin, skin factor
      elemental subroutine vertical_wll_wi_hete(wfrac, perm, length, lthfra, re, rw, skin, wi)
      real(STDD), intent(in)     :: wfrac
      real(STDD), intent(in)     :: perm
      real(STDD), intent(in)     :: length, lthfra
      real(STDD), intent(in)     :: re, rw
      real(STDD), intent(in)     :: skin
      real(STDD), intent(out)    :: wi ! well index L^3

      real(STDD), parameter      :: para = 2* ST_PI
      real(STDD) :: denominator
      
      denominator = log(re/rw) + skin
      wi = para * wfrac * perm * length * lthfra / denominator
      end subroutine vertical_wll_wi_hete
! PURPOSE:
! calculate the well-index of all the wells
! input:
! nwlls, well number
! nlys,  perforated layers number of all wells
! ilcid, perforated layers uid of all wells
! iwadj, index into perforated layers info.
! iwdir, well directions
! dwrad, well radius, ft
! dwgeo, well geofactors
! dwfra, well fractions
! dlskn, skin factor of each layer
! dlhfr, length fraction of each layer
! nclls, grid cells no.
! perm3d, =(permi, permj, permk) 3D permeability
! step3d, =(dcdx, dcdy, dcdz) grid cells step in IJK-directions
! output:
! dlwi, well index
      subroutine set_wlls_wi(nwlls, nlys, iwadj, ilcid, iwdir, dwrad, &
      dwgeo, dwfra, dlskn, dlhfr,nclls, perm3d, step3d, dlwi)
      integer, intent(in)                       :: nwlls
      integer, intent(in)                       :: nlys
      integer, dimension(nwlls+1), intent(in)   :: iwadj
      integer, dimension(nlys), intent(in)      :: ilcid
      integer, dimension(nwlls), intent(in)     :: iwdir
      real(STDD), dimension(nwlls), intent(in)  :: dwrad, dwgeo, dwfra
      real(STDD), dimension(nlys), intent(in)   :: dlskn, dlhfr
      integer, intent(in)                       :: nclls
      real(STDD), dimension(3*nclls), intent(in) :: perm3d, step3d
      real(STDD), dimension(nlys), intent(out)  :: dlwi

      integer :: widx, nlays, ly0, ly1

      do widx = 1 , nwlls
        ly0 = iwadj(widx)
        ly1 = iwadj(widx + 1) - 1
        nlays = iwadj(widx + 1) - ly0
        call set_wll_wi(nlays, ilcid(ly0:ly1), iwdir(widx), dwrad(widx),&
        dwgeo(widx), dwfra(widx), dlskn(ly0:ly1), dlhfr(ly0:ly1),       &
        nclls,perm3d, step3d, dlwi(ly0:ly1) )
      end do
      end subroutine set_wlls_wi
! calculate the well index of each well
! input:
!   nlays, perforated layers no. of a well
!   ilcid, uid of perforated layers
!   dir,    well direction
!   rw,  ft, well radius
! geofac, geofac, well in the center of the cell ~0.249
! wfrac,    well fraction                        ~1.
! dlskn, skin factor of each layer
! dlhf, perforated fraction along grid cell
! nclls, grid cells no.
! perm3d = (permi, permj, permk)
! step3d = (dcdx, dcdy, dcdz)
! output:
! wi, well index of the given well
      subroutine set_wll_wi(nlays, ilcid, dir, rw, geofac, wfrac, dlskn,&
      dlhf, nclls, perm3d, step3d, wi)
      use stmgeomet, only : get_perpend_dirs, get_geo_prop_d
      integer, intent(in)                       :: nlays
      integer, dimension(nlays), intent(in)     :: ilcid
      integer, intent(in)                       :: dir
      real(STDD), intent(in)                    :: rw
      real(STDD), intent(in)                    :: geofac
      real(STDD), intent(in)                    :: wfrac
      real(STDD), dimension(nlays), intent(in)  :: dlskn
      real(STDD), dimension(nlays), intent(in)  :: dlhf
      integer, intent(in)                       :: nclls
      real(STDD), dimension(3*nclls), intent(in):: perm3d, step3d
      real(STDD), dimension(nlays), intent(out) :: wi

      real(STDD), dimension(nlays) :: permi, permj, di, dj
      real(STDD), dimension(nlays) :: re, perm, h

      integer :: i, j
      call get_perpend_dirs(dir, i, j)
      ! block length along the well 
      call get_geo_prop_d(nclls, step3d(dir*nclls+1 : (dir+1)*nclls), nlays,&
      ilcid, h)
      ! block permeability perpendicular to the well
      call get_geo_prop_d(nclls, perm3d(i*nclls+1 : (i+1)*nclls), nlays,&
      ilcid, permi)
      call get_geo_prop_d(nclls, perm3d(j*nclls+1 : (j+1)*nclls), nlays,&
      ilcid, permj)
      ! block length perpendicular to the well
      call get_geo_prop_d(nclls, step3d(i*nclls+1 : (i+1)*nclls), nlays,&
      ilcid, di)
      call get_geo_prop_d(nclls, step3d(j*nclls+1 : (j+1)*nclls), nlays,&
      ilcid, dj)
      ! equivalent radius 
      call vertical_wll_re_hete(geofac, permi, permj, di, dj, re)
      ! permeability perpendicular to the well
      call vertical_wll_perme_hete(permi, permj, perm)
      ! well index
      call vertical_wll_wi_hete(wfrac, perm, h, dlhf, re, rw, dlskn, wi)
      end subroutine set_wll_wi
! calculate the well-index of each well at multiple layers
! If length and permeability are in ft and mD, respectively, wi in mD-ft
      subroutine wll_wi(wfrac, geofac, rw, skin, nlays,permi, permj, h, di, dj, hfrac, wi)
      real(STDD), intent(in)    :: wfrac
      real(STDD), intent(in)    :: geofac
      real(STDD), intent(in)    :: rw
      integer, intent(in)                       :: nlays
      real(STDD), dimension(nlays), intent(in)  :: skin
      real(STDD), dimension(nlays), intent(in)  :: permi, permj
      real(STDD), dimension(nlays), intent(in)  :: h, hfrac
      real(STDD), dimension(nlays), intent(in)  :: di, dj
      real(STDD), dimension(nlays), intent(out) :: wi
      
      real(STDD), dimension(nlays) :: re
      real(STDD), dimension(nlays) :: perm
      call vertical_wll_re_hete(geofac, permi, permj, di, dj, re)
      call vertical_wll_perme_hete(permi, permj, perm)
      call vertical_wll_wi_hete(wfrac, perm, h, hfrac, re, rw, skin, wi)
      wi = wi * ST_FUCONV
      call prt_wll_wi(nlays, re, perm, wi)
      end subroutine wll_wi
! PURPOSE:
! debug for wll_wi
      subroutine prt_wll_wi(n, re, perm, wi)
      integer, intent(in) :: n
      real(STDD), dimension(n), intent(in) :: re
      real(STDD), dimension(n), intent(in) :: perm
      real(STDD), dimension(n), intent(in) :: wi
      integer :: i
      write(*,'(T10, 3(1x,A10,1x))') "re(ft)    ", " perm(mD) ", &
      " wi(xxxxx)"
      do i = 1, n
        write(*, '(T10, 2(F11.3,1x), ES11.3)') re(i), perm(i), wi(i)
      enddo
      end subroutine prt_wll_wi
      end module stmwllbor
