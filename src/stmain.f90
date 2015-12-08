      program stmexe
      use stmutils 
	  use stmoutput
      implicit none

	  call st_init_memo
      write(*,'(T10,"Steam-based thermal reservoir simulator",/, &
      T14,"Jan 2016",//)') 
	  call st_reservoir_init

	  ! output into the *.out
	  call dprt_datpol
	  call st_release_memo
      end program stmexe
