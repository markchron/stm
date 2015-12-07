	  module stmutils
	  use stmheader
	  use stmdatpol
	  use stmdatnpt
      private
     
	  public :: st_init_memo, st_reservoir_init
	  contains
	  subroutine st_init_memo
      call datpol_init_scalar
	  call npt_init
	  end subroutine st_init_memo

	  subroutine st_reservoir_init
      call npt_set_properties


	  call datpol_update_index
	  end subroutine st_reservoir_init
	  end module stmutils
