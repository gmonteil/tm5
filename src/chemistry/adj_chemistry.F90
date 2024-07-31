module adj_chemistry

    use go, only : TDate

    implicit none

    public :: adj_chemistry_init, adj_chemistry_done, adj_chemistry_step
    private

    contains

        subroutine adj_chemistry_init(status)

            integer, intent(out)    :: status
            status = 10     ! Make sure it crashes, until the routine is properly implemented

        end subroutine adj_chemistry_init


        subroutine adj_chemistry_step(region, period, status)

            integer, intent(out)                    :: status
            integer, intent(in)                     :: region
            type(TDate), dimension(2), intent(in)   :: period
            status = 10     ! Make sure it crashes, until the routine is properly implemented

        end subroutine adj_chemistry_step


        subroutine adj_chemistry_done(status)

            integer, intent(out)    :: status
            status = 10     ! Make sure it crashes, until the routine is properly implemented

        end subroutine adj_chemistry_done

end module adj_chemistry