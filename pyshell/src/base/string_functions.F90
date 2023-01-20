module string_functions

implicit none

private

integer, parameter  :: long_string_len = 16384
integer, parameter  :: short_string_len = 256

public :: string_concat, string_split

contains

    function string_concat(array_of_strings, gl_str) result (long_string)

        implicit none

        ! IO
        character(len=*), dimension(:), intent(in) :: array_of_strings
        character(len=*), optional, intent(in)     :: gl_str
        character(len=long_string_len)             :: long_string

        ! LOCAL
        integer                         :: length, i, cur_index
        character(len=short_string_len) :: glue_string
        integer                         :: length_glue

        !START
        if (.not. present(gl_str) .or. len(trim(adjustl(gl_str))) == 0) then
            glue_string = ' '
            length_glue = 1
        else
            glue_string = trim(adjustl(gl_str))
            length_glue = len(trim(adjustl(glue_string)))
        end if

        long_string = ''
        cur_index = 1
        do i=1,size(array_of_strings)
            length = len(trim(adjustl(array_of_strings(i))))
            long_string(cur_index:cur_index+length-1) = trim(adjustl(array_of_strings(i)))
            cur_index = cur_index+length
            if (i < size(array_of_strings)) then
                long_string(cur_index:cur_index+length_glue-1) = glue_string(1:length_glue)
                cur_index = cur_index + length_glue
            end if
        end do

    end function string_concat

    subroutine string_split(long_string, separator, array_of_strings)

        implicit none

        ! IO
        character(len=*), intent(in)                              :: long_string
        character(len=1), intent(in)                              :: separator
        character(len=short_string_len), allocatable, intent(out) :: array_of_strings(:)

        ! LOCAL
        integer :: num_substrings, i, beg_pos, i_substring

        ! START
        num_substrings = 1
        do i = 1, len(long_string)
            if (long_string(i:i) == separator) num_substrings = num_substrings + 1
        end do

        allocate(array_of_strings(num_substrings))
        beg_pos = 1
        i_substring = 1
        do i = 1, len(long_string)
            if (long_string(i:i) == separator) then
                array_of_strings(i_substring) = trim(adjustl(long_string(beg_pos:i-1)))
                beg_pos = i+1
                i_substring = i_substring + 1
            else if (i == len(long_string)) then
                array_of_strings(i_substring) = trim(adjustl(long_string(beg_pos:i)))
                beg_pos = i
                i_substring = i_substring + 1
            end if
        end do

    end subroutine string_split

end module string_functions
