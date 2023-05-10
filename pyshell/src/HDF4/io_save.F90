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

  use io_save_hdf, only : savehdf
  use io_save_hdf, only : readhdf
  use io_save_hdf, only : readhdfmmr
  use io_save_hdf, only : readhdfmmix
  use io_save_hdf, only : readtm3hdf

  use io_save_netcdf, only : savenetcdf
  use io_save_netcdf, only : readnetcdf
  use io_save_netcdf, only : save_filename

  implicit none

  private

  public :: readnetcdf, savenetcdf, save_filename
  public :: readhdfmmr, readhdfmmix, readtm3hdf

end module io_save
