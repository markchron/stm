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
! calculate the well-index of each well at multiple layers
! If length and permeability are in ft and mD, respectively, wi in mD-ft
      subroutine wll_wi(wfrac, geofac, rw, skin, nlays,permi, permj, h, di, dj, hfrac, wi)
      real(STDD), intent(in)    :: wfrac
      real(STDD), intent(in)    :: geofac
      real(STDD), intent(in)    :: rw
      real(STDD), intent(in)    :: skin
      integer, intent(in)                       :: nlays
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
