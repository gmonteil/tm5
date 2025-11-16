#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
#include "tm5.inc"

module adj_chemistry

    use go,         only : TDate, gol, goerr
    use chemistry,  only : chemistry_step

    implicit none

    public :: adj_chemistry_init, adj_chemistry_done, adj_chemistry_step
    private

    character(len=*), parameter   :: mname = 'adj_chemistry'

    contains

        subroutine adj_chemistry_init(status)
            integer, intent(out)    :: status
            status = 0
        end subroutine adj_chemistry_init


        subroutine adj_chemistry_step(region, period, status)

            integer, intent(out)                    :: status
            integer, intent(in)                     :: region
            type(TDate), dimension(2), intent(in)   :: period
            character(len=*), parameter             :: rname = mname//'/adj_chemistry_step'

            status = 0
            call chemistry_step(region, period, status)
            IF_NOTOK_RETURN(status=1)

        end subroutine adj_chemistry_step


        subroutine adj_chemistry_done(status)
            integer, intent(out)    :: status
            status = 0
        end subroutine adj_chemistry_done

end module adj_chemistry