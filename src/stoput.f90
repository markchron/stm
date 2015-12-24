	  module stmoutput
      use stmheader
      use stmdatpol
	  ! output format of single number
	  ! 1, integer
	  ! 2, float
	  ! 3, scientific float
	  character(5), dimension(3) 		:: valfmt
	  ! output format of a line of numbers
	  character(BUF_LEN), dimension(3)  :: linefmt
	  ! number of values per line
	  integer, dimension(3)				:: novapl

      interface dprt_arr
          module procedure dprt_arr_i
          module procedure dprt_arr_d
      end interface dprt_arr
	  contains
! PURPOSE:
! create the .out file
	  subroutine dprt_create

      call dprt_create_log

      if(nprocs > 1 ) then
	    write(stofnm, '(A,"_",I4.4,A)') stifnm(1 : stiftrunc), rank, ".out"
      else
	    write(stofnm, '(A,A)') stifnm(1 : stiftrunc), ".out"
	  endif

      OPEN(UNIT = FUNIT_OUT, FILE=TRIM(stofnm), STATUS = 'REPLACE', FORM='FORMATTED', IOSTAT=stERRs(5))

      end subroutine dprt_create
! PURPOSE:
! create the .log file output unit. Only master processor has access to log file
! >> screen
! or
! >> stifnm(1:stiftrunc).log 
      subroutine dprt_create_log
      use iso_fortran_env, only : stderr => error_unit
      if(rank == MASTER) then
          FUNIT_LOG = stderr
          call st_soft_info
          call dprt_omp_system_info
      endif
      return
      end subroutine dprt_create_log
! PURPOSE:
! simply print out the array data  >> stdout
      subroutine dprt_datpol
	  character(BUF_LEN) 					:: lfmt
      call setupfmt
      call dprt_csr(Ngcll, Ngedges, ivadj, iadjncy, detrans, "transmissibility", 1)
      call dprt_csr(Nwidx, Nlyid, ilyadj, icidly, dlywitrs, "well index", 1)

	  lfmt = '(A,":",1x' // trim(valfmt(1)) //')'
	  write(FUNIT_OUT, lfmt) "Active cells no.", Nlacs
	  call dprt_arr(ictind(1:Ngcll), Ngcll, "cell type index:")
      call dprt_arr(icdist, Ngcll, "cell dist domain:")
      call dprt_arr(iwdist, Nwidx, "wellbore dist domain:")
	  call dprt_arr(dcdx, Ngcll, "cell steps along X:")
	  call dprt_arr(dcdy, Ngcll, "cell steps along Y:")
	  call dprt_arr(dcdz, Ngcll, "cell steps along Z:")
	  call dprt_arr(dchtop, Nxyplane, "cell top depth:")

	  call dprt_arr(dcpor, Ngcll, "Porosity:")
	  call dprt_arr(dcpermi, Ngcll, "Permeability I (mD)")
	  call dprt_arr(dcpermj, Ngcll, "Permeability J:")
	  call dprt_arr(dcpermk, Ngcll, "Permeability K:")
      end subroutine dprt_datpol
! set up the output format of each line	  
	  subroutine setupfmt
      integer :: i			  
	  valfmt = (/"I6   ", "F8.3 ", "G10.3"/)
      novapl = (/16, 10, 10 /)
	  do i =  1, 3
         if(novapl(i) > Nxd) novapl(i) = Nxd		  
        call dynlfmt(linefmt(i), novapl(i), valfmt(i))
      enddo
      end subroutine setupfmt
	  subroutine dynlfmt(lfmt, n, vfmt)
      character(*), intent(out) 		:: lfmt
	  integer, intent(in) 				:: n
	  character(*), intent(in) 			:: vfmt
	  write(lfmt, '( "(", I3, "(",A,",1x))" )') n, trim(vfmt)
	  end subroutine dynlfmt
      subroutine dynlfmtspcs(lfmt, n, vfmt, sps)
      character(*), intent(out) 		:: lfmt
	  integer, intent(in) 				:: n
	  character(*), intent(in) 			:: vfmt
      integer, intent(in)               :: sps
    
      character(20) ::nchar, schar
      write(nchar, '(I20)') n
      write(schar, '(I20)') sps
      lfmt ="(T"//trim(ADJUSTL(schar))//","//trim(ADJUSTL(nchar))//"("//trim(vfmt)//",1x))"
      end subroutine dynlfmtspcs
! PURPOSE:
! print out array line by line
!	  subroutine dprt_arr(arr,n,amsg, sci)
!      character(*), intent(in) 				:: amsg
!      integer, intent(in) 					:: n
!	  class(*), dimension(n), intent(in) 	:: arr
!	  integer, optional, intent(in) 		:: sci
!
!	  integer :: i, ni, npl
!	  character(BUF_LEN) 					:: lfmt
!	  lfmt = '(A,2x,"(size:"' // trim(valfmt(1)) // '" )")'
!	  write(FUNIT_OUT,lfmt) amsg, n
!	  
!	  select type (arr) 
!	  type is (real(STDD))
!		if(PRESENT(sci)) then
!			npl = novapl(3)
!			lfmt = linefmt(3)
!		else
!			npl = novapl(2)
!			lfmt = linefmt(2)
!	    endif	
!		ni = n / npl
!		do i= 1, ni
!          write(FUNIT_OUT, lfmt) arr((i-1)*npl+1 : i*npl)
!        enddo
!     	i = n - ni * npl 
!        if(i>0) write(FUNIT_OUT, lfmt) arr(ni*npl+1 : n)
!      type is (integer)
!		 npl = novapl(1)
!		 lfmt = linefmt(1)
!		 ni = n / npl
!		 do i = 1, ni
!           write(FUNIT_OUT, lfmt) arr( ((i-1)*npl+1) : (i* npl))
!         enddo
!         i = n - ni * npl
!         if ( i > 0 ) write(FUNIT_OUT, lfmt) arr(ni*npl + 1 : n)
!      end select
!      end subroutine dprt_arr
      subroutine dprt_arr_i(arr,n,amsg,sci)
      character(*), intent(in)              :: amsg
      integer, intent(in)                   :: n
      integer, dimension(:), intent(in)     :: arr
      integer, optional, intent(in)         :: sci
      integer :: i, ni, npl
      character(BUF_LEN)                    :: lfmt
      lfmt = '(A,2x,"(size:"' // trim(valfmt(1)) // '" )")'
      write(FUNIT_OUT,lfmt) amsg, n

      if(PRESENT(sci)) return

      npl = novapl(1)
      lfmt = linefmt(1)
      ni = n / npl
      do i = 1, ni
      write(FUNIT_OUT, lfmt) arr((i-1)*npl + 1 : i*npl)
      enddo
      i = n - ni * npl
      if( i> 0) write(FUNIT_OUT, lfmt) arr(ni*npl + 1 : n)
      end subroutine dprt_arr_i
      subroutine dprt_arr_d(arr, n, amsg, sci)
      character(*), intent(in)              :: amsg
      integer, intent(in)                   :: n
      real(STDD),dimension(:), intent(in)   :: arr
      integer, optional, intent(in)         :: sci
      integer :: i, ni, npl
      character(BUF_LEN)                    :: lfmt
      lfmt = '(A,2x,"(size:"' // trim(valfmt(1)) // '" )")'
      write(FUNIT_OUT,lfmt) amsg, n

      if(PRESENT(sci)) then
          npl = novapl(3)
          lfmt = linefmt(3)
      else
          npl = novapl(2)
          lfmt = linefmt(2)
      endif
      ni = n / npl
      do i = 1, ni
      write(FUNIT_OUT, lfmt) arr((i-1)*npl + 1 : i *npl)
      enddo
      i = n - ni * npl
      if(i > 0 ) write(FUNIT_OUT, lfmt) arr(ni*npl+1 : n)
      end subroutine dprt_arr_d

! print out a CSR format matrix
! input:
! sci, scientific print
      subroutine dprt_csr(n,nz,ai,aj,a, amsg, sci)
      use stmdatpol, only : maxcol_csr
      character(*), intent(in)          :: amsg
      integer, intent(in)               :: n, nz
      integer, dimension(n+1), intent(in)   :: ai
      integer, dimension(nz), intent(in)    :: aj
      class(*), dimension(nz), intent(in)   :: a
      integer, optional, intent(in)         :: sci
      ! numbers per line
      integer        ::i, npl 
      character(BUF_LEN)        :: tfmt, lfmt
      tfmt = '(A,2x,"(row:", I9,2x,"nzero:", I20, ")")'
      write(FUNIT_OUT, tfmt) amsg, n, nz
      npl = maxcol_csr(n,ai)
      call dynlfmt(lfmt, npl+1, "2xI9")
      tfmt = repeat(' ', BUF_LEN)
      select type (a)
      type is (real(STDD))
        if(PRESENT(sci)) then
            call dynlfmtspcs(tfmt, npl, "ES11.3", 16)
            do i = 1, n
            write(FUNIT_OUT, lfmt) i, aj(ai(i) : ai(i+1) -1 )
            write(FUNIT_OUT, tfmt) a( ai(i) : ai(i+1) -1 )
            enddo
        else
            call dynlfmtspcs(tfmt, npl, "G11.3", 16)
            do i = 1, n
            write(FUNIT_OUT, lfmt) i, aj(ai(i) : ai(i+1) -1 )
            write(FUNIT_OUT, tfmt) a( ai(i) : ai(i+1) -1 )
            enddo
         endif
        type is(integer)
            call dynlfmtspcs(tfmt, npl, "I11  ", 16)
            do i = 1, n
            write(FUNIT_OUT, lfmt) i, aj(ai(i) : ai(i+1) -1 )
            write(FUNIT_OUT, tfmt) a( ai(i) : ai(i+1) -1 )
            enddo
        end select
      end subroutine dprt_csr
! PURPOSE:
! print out the omp info
      subroutine dprt_omp_system_info
      use omp_lib
      write(FUNIT_LOG, 6800) omp_get_num_procs(), omp_get_max_threads(), omp_get_num_threads()
 6800 format(1x, 55('-'), /, T21, 'omp_get_num_procs', T45, I10,/, &
      T21, 'omp_get_max_threads', T45, I10,/,  &
      T21, 'omp_get_num_threads', T45, I10,/, 1x, 55('-') )
      end subroutine dprt_omp_system_info
! PURPOSE:
! print software info
     subroutine st_soft_info
      write(FUNIT_LOG,'(T10,"Steam-based thermal reservoir simulator",/, &
      T30,"Jan 2016",//)') 
     end subroutine st_soft_info
! PURPOSE
! check the stErrs info.
	  subroutine dprt_errmsg
      integer :: i
	  i = sum(stErrs)
	  if(i /= 0 ) then 
		  write(FUNIT_LOG, '(TR10, "There are errors in stErrs." )')
	  else
		  write(FUNIT_LOG, '("STM")')
	  endif
	  end subroutine dprt_errmsg

      end module stmoutput
