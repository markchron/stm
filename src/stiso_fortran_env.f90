      module iso_fortran_env
      ! noninstrinsic version for Fortran of Linux
      ! see subclause 13.8.2 of the fortran 2003 standard.

      implicit NONE
      public
      integer, parameter :: Character_Storage_Size = 8
      integer, parameter :: error_unit             = 0
      integer, parameter :: File_Storage_Size      = 8
      integer, parameter :: input_unit             = 5
      integer, parameter :: iostat_end             = -1
      integer, parameter :: iostat_eor             = -2
      integer, parameter :: numeric_storage_size   = 32
      integer, parameter :: output_unit            = 6

      end module iso_fortran_env

