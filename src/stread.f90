	  module stmdatnpt
      use stmheader
	  use stmdatpol
      integer, parameter 		:: I_UNIT_5 = 15
	  character(4), parameter 	:: I_DAT_SUFFIX='.dat'
	  contains
! pre-read the data size | structure related info.
! predefine all these  values by revise the program
! define the global variables
      subroutine npt_init
      Nxd = 9
	  Nyd = 9
      Nzd = 4
	  Nwidx = 3
	  NLYID = 9  ! defined after go through the datset, here, 3 vertical wells
	  fnstec = 0 		! pre-read:  discretized format, 5 or 9 points
	  metric = 0 		! field units

      Ncidx = 2 ! components no.
      Npidx = 3 ! phases no.
      Nevents = 5
	  call npt_update
      end subroutine npt_init
! update the data size related info 
	  subroutine npt_update
	  Nxyplane = Nxd * Nyd
      Ngcll = Nxyplane * Nzd ! total grid cells
      Ngidx = Ngcll + Nwidx  ! total blocks

	  Neqn = Ncidx + 1 + Npidx + 1
	  Npeqn = Ncidx + 1
	  end subroutine npt_update

	  subroutine npt_set_properties
      integer :: i, k
      ictind = (/ ( 5, 2, 2, 2, 2, 2, 2, 2, 5,  &
			  0, 4, 1, 1, 1, 1, 1, 3, 0, &
			  0, 0, 4, 1, 1, 1, 3, 0, 0,  &
			  0, 0, 0, 4, 1, 3, 0, 0, 0, &
			  0, 0, 0, 0, 6, 0, 0, 0, 0, &
			  (0, i=1, 36), k=1, Nzd) /)
	  dcdx = 29.17 ! ft
	  dcdy = 29.17 ! ft
	  dcdz = (/ (10., i=1, Nxyplane),  (20., i=1, Nxyplane), &
			  (25., i=1, Nxyplane),  (25., i=1, Nxyplane)  /) ! ft, KVAR
	  dchtop = 1500 ! ft
	  
	  dcpor = 0.3
	  dcpermi  =  (/ (2000., i=1, Nxyplane),  (500., i=1, Nxyplane), &
			  (1000., i=1, Nxyplane),  (2000., i=1, Nxyplane)  /) ! Darcy, KVAR
	  dcpermj = dcpermi
	  dcpermk = dcpermi * 0.5

      !  UID   ff    skin  status       connection    
      ! 1 1 4   1.     0.     open        flow-from  'surface'
      !
      !  UID   ff    skin  status       connection    
      ! 9 1 1   1.     0.     open        flow-to  'surface'
      ! 9 1 2   1.     0.     open        flow-to  1
      ! 9 1 3   1.     0.     open        flow-to  2      
      ! 9 1 4   1.     0.     open        flow-to  3
      !
      !  UID   ff    skin  status       connection    
      ! 5 5 1   1.     0.     open        flow-to  'surface'
      ! 5 5 2   1.     0.     open        flow-to  1
      ! 5 5 3   1.     0.     open        flow-to  2      
      ! 5 5 4   1.     0.     open        flow-to  3
      ixly = (/ 1, (9, i=1, 4), (5, i=1, 4) /)
      iyly = (/ 1, (1, i=1, 4), (5, i=1, 4) /)
      izly = (/ 4, (i, i=1, 4), (i, i=1, 4) /)
      ialy = (/ 0, (i, i=0, 3), (i, i=0, 3) /)
      ! The perforation info of (widx)-well is stored in ixly(iyly,izly,...)
      ! starting at index ilyadj(widx) and ending at (but not including) index
      ! ilyadj(widx+1) (i.e. ixly(ilyadj(widx) : ilyadj(widx+1)-1) store all the
      ! perforated UI index of well 'widx'
      ilyadj = (/ 1, 2, 6, 10 /) 
      iwdir = (/ (2, i=1, Nwidx) /) ! vertical well
      dwllrad = (/ (0.3, i=1, Nwidx) /) ! ft
      dwllgeof = (/ (0.249d0, i=1, Nwidx) /)
      dwllfrac = (/ (1.d0, i=1, Nwidx) /)

      dlylenfr = (/ (1.d0, i=1, NLYID) /)
      dlyskin = (/ (0.d0, i=1, NLYID) /)
     
      ! data|events offset with initial, days
      devnts =(/ (5.*i, i = 0, Nevents-1) /)
	  end subroutine npt_set_properties
	  
! -------------------------------------------------------------------------
! PURPOSE
! read or get the input file name by screen input or after command line
! for example
! ./exe test.dat --> infile='test.dat'//repeat(' ',xxx)
!                    trunc= 4, which represents the last character
!                    without the suffix. 
! if need infile length without space characters, LEN_TRIM(infile)
! It also broadcast the input file name and its position without suffix
! to the other processors, if it called in a MPI environment.
! MPI: each processor should check the input file existence. 
      subroutine get_ifile
      if(rank .eq. MASTER) then
          call serial_get_ifile
      endif
      end subroutine get_ifile
! PURPOSE:
! check the command line with input file or not. If no input file, allow
! an file input from keyboard.
! DISABLED FUNCTION: interpret the input line, by file name and stored
! path, when the commented subroutine 'serial_intepret_string is
! actived. 
      subroutine serial_get_ifile
      character(len=BUF_LEN) buf
      if(IARGC() .eq. 0) then
        WRITE(*, 6801)
        READ(*, 6802) buf
      elseif(IARGC() .ge. 1) then
        call GETARG(1, buf)
      endif
      stifnm = TRIM(buf)
      stiflen = LEN_TRIM(buf)

      call serial_truncate_filename
 6801 FORMAT(T5,'Please enter the data file name :')
 6802 FORMAT(A)
      end subroutine serial_get_ifile
! PURPOSE:
! Truncate the file without input data file suffix. For example:
! (1) infile = test.dat,  --> trunc = 4
! (2) infile = filename   --> trunc = 8
! (3) infile=test.da      --> trunc = 7
      subroutine serial_truncate_filename
      stiftrunc = SCAN( stifnm(1 : stiflen), '.', .TRUE.)
      if(stiftrunc /= 0) then
        if( stifnm(stiftrunc : stiflen) == I_DAT_SUFFIX)then 
          stiftrunc = stiftrunc - 1
          return
        endif
      endif 
! stiftrunc == 0, can not find '.' OR stiftrunc/=0 (find '.') but '.xxx'
! /='.dat'
      stiftrunc= stiflen
      return
      end subroutine serial_truncate_filename
! PURPOSE:
! check the input file existence.      
      subroutine serial_verify_ifile
      OPEN(UNIT= I_UNIT_5, FILE=stifnm(1:stiflen), STATUS='OLD', ACTION='READ', IOSTAT=stErrs(3))
      if(stErrs(3) /= 0) then! input file does not exist
		  stErrs(4) = 1 ! input file does not exist
          WRITE(*,6802) rank, stifnm
      else ! input file exists
          WRITE(FUNIT_OUT, 6804) rank, stifnm(1 : stiflen)
      endif! file exits?

 6802 FORMAT(/,T10,'Error 0004: reported from processor ',I4,/,T20,'File "',A,'" does not exist.')
 6804 FORMAT(T18,'CHECK:',2x,'Input file@',I4,':',T45, A)
      end subroutine serial_verify_ifile
!PURPOSE:
! intepret the string. Truncate the PATH of the given file. 
! For example, input file name: string= /home/xxx/workspace/test/test.dat
! after this subroutine, return, pos = 25
!                               infile = test.dat
!                               lth_ifile= 8 = 33 - 25
      subroutine serial_intepret_string(string, pos, infile,lth_ifile)
      character(*), intent(in) :: string
      integer, intent(out) :: pos
      character(*), intent(out) :: infile
      integer, intent(out) :: lth_ifile
      integer lth
      
      lth = LEN_TRIM(string)
      pos = SCAN(string,'\/',.TRUE.)
      if(pos /= 0) then
        lth_ifile = lth - pos
      else ! pos == 0
        lth_ifile = lth
      endif
      if(lth_ifile .gt. BUF_LEN) then
        stErrs(3) = 1
        WRITE(*, 6803) string(pos+1 : lth)
      else
        infile = string(pos+1:lth)
      endif
      
      if(pos /=0) WRITE(FUNIT_OUT, 6804) string(1:pos)
      
 6803 FORMAT(/,T10,'Error:',1x,'Input file',/, T3,'"',A,'"',/, T21,'is too long, and its name is truncated.')
 6804 FORMAT('0','Input file locates at ',2x,A)
      end subroutine serial_intepret_string	  
	  end module stmdatnpt
