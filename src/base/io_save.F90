!###############################################################################
!
! contains routines to read emissions and
! to read and write the main model state from/to file
!
!### macro's ###################################################################
!
#include "tm5.inc"
!
!###############################################################################

module io_save

  use io_save_netcdf, only : savenetcdf
  use io_save_netcdf, only : readnetcdf
  use io_save_netcdf, only : save_filename

  implicit none

  private

  public :: readnetcdf, savenetcdf, save_filename

end module io_save
