c
c Module to include "ununits.inc" provided with the UDUnits package.
c Thanks to its Fortran-77 style it could not be included in the
c 'udunits.F90' module directly.
c
c Most important part of the include file is the definition of 'UD_POINTER',
c which specfies how in Fortran a pointer to a uni-data type
c should be defined. Depending on the installed 'udunits' library,
c this is either 'integer*4' or 'integer*8' :
c
c       #define UD_POINTER integer*8
c
c With this macro the creation function is defined:
c
c       UD_POINTER utmake
c

      module udunits_inc
      
        public
      
#include "udunits.inc"

        ! set the kind value given the kind of the UTMake function:
        integer, parameter  ::  UD_POINTER_KIND = kind(utmake)

      end module udunits_inc

