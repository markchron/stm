      program stmexe
      use stmutils 
	  use stmoutput
      implicit none

	  call st_init_memo
	  call st_reservoir_init
      call st_dist
	  ! output into the *.out
	  call dprt_datpol
	  call st_release_memo
      end program stmexe
