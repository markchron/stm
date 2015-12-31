      module stmvtkxml
      use stmheader, only : STDD, BUF_LEN
      contains
! PURPOSE:
! print software info
      subroutine prt_soft_info(funit)
      integer, intent(in)           :: funit
      write(funit,'(T10,"Steam-based thermal reservoir simulator",/, &
      T30,"Jan 2016",//)') 
      end subroutine prt_soft_info
! PURPOSE:
! print out the omp info
      subroutine prt_omp_system_info(funit)
      use omp_lib
      integer, intent(in)           :: funit
      write(funit, 6800) omp_get_num_procs(), omp_get_max_threads(), omp_get_num_threads()
 6800 format(1x, 55('-'), /, T21, 'omp_get_num_procs', T45, I10,/, &
      T21, 'omp_get_max_threads', T45, I10,/,  &
      T21, 'omp_get_num_threads', T45, I10,/, 1x, 55('-') )
      end subroutine prt_omp_system_info
! PURPOSE:
! set up dynamic line format string
	  subroutine dynlfmt(lfmt, n, vfmt)
      character(*), intent(out) 		:: lfmt
	  integer, intent(in) 				:: n
	  character(*), intent(in) 			:: vfmt
	  write(lfmt, '( "(", I3, "(",A,",1x))" )') n, trim(vfmt)
	  end subroutine dynlfmt
! set up dynamic line format string with blank ahead
! sps, specify the blank length
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
! set up dynamic line format with title name ahead
! tfmt, specifies the title string fmt
      subroutine dynlfmttit(lfmt, tfmt, n, vfmt)
      character(*), intent(out)             :: lfmt
      character(*), intent(in)              :: tfmt
      integer, intent(in)                   :: n
      character(*), intent(in)              :: vfmt
      
      character(20) :: nchar
      write(nchar, '(I20)') n
      lfmt = "("//trim(tfmt)//","//trim(ADJUSTL(nchar))//"("//trim(vfmt)//",1x))"
      end subroutine dynlfmttit
! PURPOSE:
! print out a CSR format matrix
      subroutine prt_mtx_csr_d(funit, n, nnz, ai, aj, a, vfmt)
      integer, intent(in)               :: funit
      integer, intent(in)               :: n, nnz
      integer, dimension(n+1), intent(in)       :: ai
      integer, dimension(nnz), intent(in)       :: aj
      real(STDD), dimension(nnz), intent(in)      :: a
      character(*), intent(in)                  :: vfmt

      integer       :: npl, i
      character(BUF_LEN)        :: tfmt, lfmt, nchar
      npl = MAXVAL( ai(2:n+1) - ai(1:n) )
      write(nchar, '(I20)') npl
      nchar = ADJUSTL(nchar)
      lfmt ="(I9,1x,"//trim(nchar)//"(I9,2x))"
      tfmt ="(T10,"  //trim(nchar)//"("//trim(vfmt)//",1x))"
      do i=1, n
       write(funit, lfmt) i, aj(ai(i) : ai(i+1) -1)
       write(funit, tfmt)    a (ai(i) : ai(i+1) -1)
      enddo
      end subroutine prt_mtx_csr_d
      subroutine prt_mtx_csr_i(funit, n, nnz, ai, aj, a, vfmt)
      integer, intent(in)               :: funit
      integer, intent(in)               :: n, nnz
      integer, dimension(n+1), intent(in)       :: ai
      integer, dimension(nnz), intent(in)       :: aj
      integer, dimension(nnz), intent(in)      :: a
      character(*), intent(in)                  :: vfmt

      integer       :: npl, i
      character(BUF_LEN)        :: tfmt, lfmt, nchar
      npl = MAXVAL( ai(2:n+1) - ai(1:n) )
      write(nchar, '(I20)') npl
      nchar = ADJUSTL(nchar)
      lfmt ="(I9,1x,"//trim(nchar)//"(I9,2x))"
      tfmt ="(T10,"  //trim(nchar)//"("//trim(vfmt)//",1x))"
      
      do i=1, n
       write(funit, lfmt) i, aj(ai(i) : ai(i+1) -1)
       write(funit, tfmt)    a (ai(i) : ai(i+1) -1)
      enddo
      end subroutine prt_mtx_csr_i
      subroutine prt_csr(funit, n, nnz, ai, aj, a, amsg, sci)
      integer, intent(in)               :: funit
      integer, intent(in)               :: n, nnz
      integer, dimension(n+1), intent(in)       :: ai
      integer, dimension(nnz), intent(in)       :: aj
      class(*), dimension(nnz), intent(in)      :: a
      character(*), intent(in)          :: amsg
      integer, intent(in)               :: sci

      character(BUF_LEN)            :: tfmt
      tfmt = '(A,2x,"(row:", I9, 2x,"nzero:", I20 ")")'
      write(funit, tfmt) amsg, n, nnz
      select type(a)
      type is (real(STDD)) 
          if(sci == 1) then
              call prt_mtx_csr_d(funit, n, nnz, ai, aj, a, "ES11.3")
          else
              call prt_mtx_csr_d(funit, n, nnz, ai, aj, a, "G11.3")
          endif
      type is (integer)
          call prt_mtx_csr_i(funit, n, nnz, ai, aj, a, "I11")
      end select
      end subroutine prt_csr
! PURPOSE:
! print out array line by line
      subroutine prt_arr_d(funit, arr,n, lfmt, npl)
      integer, intent(in)                   :: funit
      integer, intent(in)                   :: n
      real(STDD),dimension(:), intent(in)   :: arr
      character(*), intent(in)              :: lfmt
      integer, intent(in)                   :: npl

      integer :: i, ni
      ni = n / npl
      do i = 1, ni
      write(funit, lfmt) arr( (i-1)*npl + 1: i*npl )
      enddo
      i = n - ni * npl
      if( i>0 ) write(funit, lfmt) arr( ni*npl+1 : n)
      end subroutine prt_arr_d
      subroutine prt_arr_i(funit, arr, n, lfmt, npl)
      integer, intent(in)                   :: funit
      integer, intent(in)                   :: n
      integer,dimension(:), intent(in)      :: arr
      character(*), intent(in)              :: lfmt
      integer, intent(in)                   :: npl

      integer :: i, ni
      ni = n / npl
      do i = 1, ni
      write(funit, lfmt) arr( (i-1)*npl + 1: i*npl )
      enddo
      i = n - ni * npl
      if( i>0 ) write(funit, lfmt) arr( ni*npl+1 : n)
      end subroutine prt_arr_i
!	  subroutine prt_arr(funit, arr,n, lfmt, npl)
!      integer, intent(in)                  :: funit
!      integer, intent(in) 					:: n
!	  class(*), dimension(n), intent(in) 	:: arr
!      character(*), intent(in) 			:: lfmt
!	  integer, intent(in) 		            :: npl
!
!	  integer :: i, ni
!	  
!	  ni = n / npl
!	  select type (arr) 
!	  type is (real(STDD))
!		do i= 1, ni
!          write(funit, lfmt) arr((i-1)*npl+1 : i*npl)
!        enddo
!     	i = n - ni * npl 
!        if(i>0) write(funit, lfmt) arr(ni*npl+1 : n)
!      type is (integer)
!		 do i = 1, ni
!           write(FUNIT, lfmt) arr( ((i-1)*npl+1) : (i* npl))
!         enddo
!         i = n - ni * npl
!         if ( i > 0 ) write(FUNIT, lfmt) arr(ni*npl + 1 : n)
!      end select
!      end subroutine prt_arr
      subroutine prt_scalar(funit, msg, var)
      integer, intent(in)            :: funit
      character(*), intent(in)       :: msg
      class(*), intent(in)           :: var
      select type (var)
      type is (integer) 
         write(funit, '(T2,A,":",T42,I6)') msg, var
      type is(real(STDD))
         write(funit, '(T2,A,":",T42,ES9.2)') msg, var
      end select
      end subroutine prt_scalar
! print out matrix with column and row titles
! Note: the matrix is printed out column by column since Fortran stores matrix
! in column-major
      subroutine prt_mtx_wtitles_i(funit, amsg, m,n, mtx, colt, rowt)
      integer, intent(in)                   :: funit
      character(*), intent(in)              :: amsg
      integer, intent(in)                   :: m,n
      integer, dimension(m, n), intent(in)  :: mtx
      character(*), dimension(n), intent(in) :: colt
      character(*), dimension(m), intent(in) :: rowt
      integer       :: c 
      character(BUF_LEN) lfmt
      write(funit, '(A)') amsg
      call dynlfmtspcs(lfmt, m, 'A12', 12)
      write(funit, lfmt) rowt
      call dynlfmttit(lfmt, 'A10,2x', n, 'I12')
      do c = 1, n
        write(funit, lfmt) colt(c), mtx(:, c)
      enddo

      end subroutine prt_mtx_wtitles_i
      end module stmvtkxml
