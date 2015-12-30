	  module stmoutput
      use stmheader
      use stmdatpol
      use stmvtkxml
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
! PURPOSE:
! print out array line by line
      subroutine dprt_arr_i(arr,n,amsg,sci)
      character(*), intent(in)              :: amsg
      integer, intent(in)                   :: n
      integer, dimension(:), intent(in)     :: arr
      integer, optional, intent(in)         :: sci
      character(BUF_LEN)                    :: lfmt
      lfmt = '(A,2x,"(size:"' // trim(valfmt(1)) // '" )")'
      write(FUNIT_OUT,lfmt) amsg, n

      if(PRESENT(sci)) return

      call prt_arr_i(FUNIT_OUT, arr, n, linefmt(1), novapl(1))
      end subroutine dprt_arr_i
      subroutine dprt_arr_d(arr, n, amsg, sci)
      character(*), intent(in)              :: amsg
      integer, intent(in)                   :: n
      real(STDD),dimension(:), intent(in)   :: arr
      integer, optional, intent(in)         :: sci
      character(BUF_LEN)                    :: lfmt
      lfmt = '(A,2x,"(size:"' // trim(valfmt(1)) // '" )")'
      write(FUNIT_OUT,lfmt) amsg, n

      if(PRESENT(sci)) then
          call prt_arr_d(FUNIT_OUT, arr, n, linefmt(3), novapl(3))
      else
          call prt_arr_d(FUNIT_OUT, arr, n, linefmt(2), novapl(2))
      endif
      end subroutine dprt_arr_d

! print out a CSR format matrix
! input:
! sci, scientific print
      subroutine dprt_csr(n,nz,ai,aj,a, amsg, sci)
      character(*), intent(in)          :: amsg
      integer, intent(in)               :: n, nz
      integer, dimension(n+1), intent(in)   :: ai
      integer, dimension(nz), intent(in)    :: aj
      class(*), dimension(nz), intent(in)   :: a
      integer, optional, intent(in)         :: sci
      
      integer :: isc
      if(PRESENT(sci)) isc = 1
      call prt_csr(FUNIT_OUT, n, nz, ai, aj, a, amsg, isc)
      end subroutine dprt_csr
! PURPOSE:
! print out the omp info
      subroutine dprt_omp_system_info
      call prt_omp_system_info(FUNIT_LOG)
      end subroutine dprt_omp_system_info
! PURPOSE:
! print software info
     subroutine st_soft_info
     call prt_soft_info(FUNIT_LOG)
     end subroutine st_soft_info
! PURPOSE
! check the stErrs info.
	  subroutine dprt_errmsg
      integer :: i
      if(rank /= MASTER) return 
	  i = sum(stErrs)
	  if(i /= 0 ) then 
		  write(FUNIT_LOG, '(TR10, "There are errors in stErrs." )')
	  else
		  write(FUNIT_LOG, '("STM")')
	  endif
	  end subroutine dprt_errmsg

      end module stmoutput
