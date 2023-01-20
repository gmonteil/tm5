!
!ProTeX: 1.14-AJS
!
!BOI
!
! !TITLE:        GO  -  General Objects library
! !AUTHORS:      Arjo Segers
! !AFFILIATION:  KNMI
! !DATE:         \today
!
! !INTRODUCTION: Introduction
!
! The GO library provides general purpose entities for Fortran90 programs.
!
! In summary, the GO library provides:
!
! \begin{itemize}
!   \item
!     facilities to write (error) messages;
!   \item
!     system constants and interface to system functions;
!   \item
!     character functions;
!   \item
!     date and time manipulation;
!   \item
!     file handling without specification of unit numbers;
!   \item
!     reading from text file with comment lines;
!   \item
!     reading settings from a resource file;
! \end{itemize}
!
! The library consists of a single module to provide access to routines,
! functions, and data types.
! In addition, a few shell scripts are provided for specialized tascs.
!
!
! !INTRODUCTION: About status argument
!
!   An integer status argument is present in most routines:
!   \bv
!     subroutine aaa( status )
!
!      integer, intent(out)  ::  status
!
!      ! wrong!
!      if ( 1 /= 0 ) then
!        write (*,'("ERROR - this is a serious error ...")')
!        write (*,'("ERROR in aaa"); status=1; return
!      end if
!
!      ! ok
!      status = 0
!
!    end subroutine aaa
!  \ev
!  The value should be checked on return to test on succesful completion:
!  \bv
!    call aaa( status )
!    if (status/=0) then; write (*,'("ERROR in prog")'); status=1; return; end if
!  \ev
!  To make the make the code more readible, define a preprossor macro
!  for this check-and-trace line:
!  \bv
!    #define IF_ERROR_THEN(action) if (status/=0) then; write (*,'("ERROR in prog")'); action; return; end if
!
!    call aaa( status )
!    IF_ERROR_RETURN( status=1 )
!  \ev
!
!
! !INTRODUCTION: GO
!
!   The GO-library takes the form of a module:
!   \bv
!     use GO
!   \ev
!
!   The GO module (source file \texttt{go.f90}) is in fact a shell around a bunch of
!   sub modules.
!   The sub modules should not be used directly, but accessed through
!   the main module only.
!
!
!
! !INTRODUCTION: GO_FU
!
!   \subsection{Description}
!
!     Standard file units, might be compiler depended.
!
!     Since the standard file units might be used by other modules
!     to write (error) messages, this is the most basic of all GO modules.
!
!   \subsection{Usage}
!
!     Access the entities of this module via the main module (prefered)
!     or directly:
!     \bv
!       use GO, only : ...
!       use GO_FU, only : ...
!     \ev
!
!
!   \subsection{Constants}
!
!     \begin{itemize}
!       \item
!         Standard input, output, and error file:
!         \bv
!           integer, parameter  ::  goStdIn  = 0
!           integer, parameter  ::  goStdOut = 5
!           integer, parameter  ::  goStdErr = 6
!         \ev
!       \item
!         Range of file units that may be used by this program:
!         \bv
!           integer, parameter    ::  goFuRange(2) = (/0,999/)
!         \ev
!         Free file units selected by 'goGetFu' in module 'GO_File'
!         are within this range; if all are in use, an error is issued.
!         This range might become a variable in future, such that free
!         file units are selected in a range that is not used by other
!         parts of the code.
!     \end{itemize}
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_fu.F90} .
!
!
!
! !INTRODUCTION: GO_Print
!
!   \subsection{Description}
!
!     Tools for (error) messages.
!
!     Basic idea is to write messages to the character buffer 'gol'
!     and let the the actual printing be managed by a special routine.
!     This is especially usefull if the program runs on multiple processors:
!     \bv
!       write (gol,'("This is processor : ",i2)') myid; call goPr
!     \ev
!     With the appropriate settings, the output might be:
!     \bv
!       [00] This is processor :  0
!       [03] This is processor :  3
!       [02] This is processor :  2
!       [01] This is processor :  1
!     \ev
!
!   \subsection{Usage}
!
!     Access the entities of this module via the main module (prefered)
!     or directly:
!     \bv
!       use GO, only : ...
!       use GO_Print, only : ...
!     \ev
!
!   \subsection{Variables}
!
!     \begin{itemize}
!       \item
!         Buffer for standard output:
!         \bv
!           character(len=LONG_STR_LEN)   ::  gol
!         \ev
!     \end{itemize}
!
!   \subsection{Procedures}
!
!     \begin{itemize}
!       \item
!         Initialize printing of messages:
!         \bv
!           call GO_Print_Init( status, apply=.true., &
!                                 prompt_pe=.true., pe=0, &
!                                 trace=.false. )
!         \ev
!         Arguments:
!         \begin{itemize}
!           \item \texttt{integer, intent(inout)  ::  status}\\
!             Should be zero on input, unequal to zero in case of errors.
!           \item \texttt{logical, intent(in), optional :: apply}\\
!             If set to true (on this processor), messages are printed.
!             Use this flag to limit printing to root processor only for example.
!           \item \texttt{logical, intent(in), optional  ::  prompt_pe}\\
!             If true, the processor number is printed at the begin of each line.
!           \item \texttt{integer, intent(in), optional  ::  pe}\\
!             Processor number used if \texttt{prompt_pe} is true.
!           \item \texttt{logical, intent(in), optional  ::  trace}\\
!             Print all labels (see routine \texttt{goLabel}).
!         \end{itemize}
!       \item
!         Donialize printing:
!         \bv
!           call GO_Print_Done( status )
!         \ev
!         Arguments:
!         \begin{itemize}
!           \item \texttt{integer, intent(inout)  ::  status}\\
!             Should be zero on input, unequal to zero in case of errors.
!         \end{itemize}
!       \item
!         Print the text written to '\texttt{gol}', including prompts etc:
!         \bv
!           write (gol,'("Hello!")'); call goPr
!         \ev
!       \item
!         Write error messages:
!         \bv
!           write (gol,'("this is a serious error ...")'); call goErr
!           write (gol,'("in program")'); call goErr; status=1; return
!         \ev
!         This will produce the following output:
!         \bv
!           [00] ERROR - this is a serious error ...
!           [00] ERROR - in program
!         \ev
!         At the moment, error messages are treated the same as normal messages.
!       \item
!         Set and unset routine labels:
!         \bv
!           call goLabel( 'mymod/myroutine' )
!           ...
!           call goLabel()
!         \ev
!         The labels can be used to show in which routine info is printed
!         or where exactelly an error occured.
!         If option '\texttt{trace}' in the initialization is true,
!         labels are always printed; nice for debugging.
!         By the first command, a new label is pushed on a stack to let
!         the program known that a new part of the code is reached;
!         the second command (without arguments) indicates that the end
!         of this part is reached and the label is popped from the stack.
!         Messages printed in between the above commands will be indented;
!         the height of the stack determines the indention.\\
!         A program with three printing levels, and '\texttt{trace}' and
!         '\texttt{prompe_pe}' set to true, might produce the following output:
!         \bv
!           [00] <myprog>
!           [00]   This is the main program.
!           [00]   <mymod/aaa>
!           [00]     This is routine aaa.
!           [00]     <mymod/bbb>
!           [00]       This is routine bbb.
!           [00]     (mymod/bbb)
!           [00]   (mymod/aaa)
!           [00]   Normal end.
!           [00] (myprog)
!         \ev
!         The current label is also unset by a second call to '\texttt{goErr}'
!         if '\texttt{gol}' has not been filled in between:
!         \bv
!           subroutine aaa( status )
!
!             integer, intent(out) ::  status
!
!             call goLabel( 'mymod/aaa' )
!
!             call bbb( status )
!             if (status/=0) then; call goErr; status=1; return; end if
!
!             status=0; call goLabel()
!
!           end subroutine aaa
!
!           subroutine bbb( status )
!
!             integer, intent(out) ::  status
!
!             call goLabel( 'mymod/bbb' )
!
!             write (gol,'("this is a serious error ...")'); call goErr
!             call goErr; status=1; return
!
!             status=0; call goLabel()
!
!           end subroutine bbb
!         \ev
!         This will produce the following output (no trace):
!         \bv
!           [00]       ERROR - this is a serious error ...
!           [00]       ERROR - in mymod/bbb
!           [00]     ERROR - in mymod/aaa
!           [00]   ERROR - in myprog
!         \ev
!     \end{itemize}
!
!   \subsection{To be done}
!
!     The following options should be implemented soon, since they have
!     been implemented in older versions:
!     \begin{itemize}
!       \item option to write output (per processor?) to a file
!         instead of standard output;
!       \item error messages should be written to standar error or an error file;
!       \item there should be a simple swith to turn on message printing
!         for a while
!       \item the time spent on a labelled code could be counted;
!         a very simple but effective way to profile program performance.
!     \end{itemize}
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_print.f90} .
!
!
!
! !INTRODUCTION: GO_String
!
!   General objects for character strings.
!
!   \subsection{Usage}
!
!     Access the entities of this module via the main module (prefered)
!     or directly:
!     \bv
!       use GO, only : ...
!       use GO_String, only : ...
!     \ev
!
!   \subsection{Procedures}
!
!     \begin{itemize}
!       \item
!         Split a string into two parts:
!         \bv
!           call goSplitLine( 'ab#cd', s1, '#', s2, status )
!         \ev
!         The input string is split at the first occurance of \texttt{'#'};
!         the leading part is returned in \texttt{s1}, and the remainder
!         without \texttt{'#'} in \texttt{s2}.
!         One or both of \texttt{s1} and \texttt{s2} might be empty. \newline
!         Interface:
!         \bv
!           subroutine goSplitLine( line, s1, c, s2 )
!             character(len=*), intent(in)      ::  line
!             character(len=*), intent(out)     ::  s1
!             character(len=1), intent(in)      ::  c
!             character(len=*), intent(out)     ::  s2
!             integer, intent(inout)            ::  status
!           end subroutine goSplitLine
!         \ev
!       \item
!         Remove leading part of a string and store contents in a variable:
!         \bv
!           line = 'abc, 123'
!           call goReadFromLine( line, x, status [,sep=','] )
!         \ev
!         The character string \texttt{line} is split at the first komma
!         (or the character specified by the optional argument \texttt{sep}),
!         the contents of the leading part is read into the
!         variable \texttt{x} and the remainder is returned in \texttt{line}.
!         Currently only variables of standard type could be read. \newline
!         Interface:
!         \bv
!           subroutine goReadFromLine( line, x, status, sep )
!             character(len=*), intent(inout)              ::  line
!             <xtype>, intent(out)                         ::  x
!             integer, intent(inout)                       ::  status
!             character(len=*), intent(in), optional       ::  sep
!           end subroutine goReadFromLine
!
!           <xtype> = integer | real | logical | character(len=*)
!         \ev
!       \item
!         Reads the value assigned to a name from a character line:
!         \bv
!           bb = 'default'
!           call goVarValue( 'aa=1;bb=xyz;cc=U123', ';', 'bb', '=', bb, status )
!         \ev
!         The character line is split at the specified seperations (here ';'),
!         and examined on the occurance of the specified name (here 'bb')
!         followed by an assignment character ('=') .
!         The value is storred in the last variable before the status argument.
!         If the name is not found, the value remains what it was.
!         Interface:
!         \bv
!           subroutine goVarValue( line, sep, var, is, val, status )
!             character(len=*), intent(in)     ::  line
!             character(len=1), intent(in)     ::  sep
!             character(len=*), intent(in)     ::  var
!             character(len=1), intent(in)     ::  is
!             <xtype>, intent(inout)           ::  x
!             integer, intent(inout)           ::  status
!           end subroutine goVarValue( line, sep, var, is, val, status )
!
!           <xtype> = integer | character(len=*)
!
!           Return status:
!             <0  :  variable not found, val remains the same
!              0  :  variable found, val reset
!             >0  :  error
!         \ev
!       \item
!         Convert integer value to a character string:
!         \bv
!           s = goNum2Str( i [,fmt='(i6)'] )
!         \ev
!         Returns a length 6-character string with the representation of the
!         integer value \texttt{i} in the first characters.
!         An optional format could be provided, for example to include
!         leading zero's.
!         Interface:
!         \bv
!           character(len=6) function goNum2Str( i, fmt )
!             integer, intent(in)                     ::  i
!             character(len=*), intent(in), optional  ::  fmt
!           end function goNum2Str
!         \ev
!       \item
!         Convert to upper/lower case:
!         \bv
!           s2 = goUpCase( s1 )
!           s2 = goLoCase( s1 )
!         \ev
!         Interfaces:
!         \bv
!           function goUpCase( s )
!             character(len=*), intent(in)    ::  s
!             character(len=len(s))           ::  goUpCase
!           end function goUpCase
!
!           function goLoCase( s )
!             character(len=*), intent(in)    ::  s
!             character(len=len(s))           ::  goLoCase
!           end function goLoCase
!         \ev
!       \item
!         Write character keys formed by text and number:
!         \bv
!           call goWriteKeyNum( key, 'sh', 159 )   !  sh159
!         \ev
!         Interface:
!         \bv
!           subroutine WriteKeyNum( res, key, num )
!             character(len=*), intent(out)    ::  res
!             character(len=*), intent(in)     ::  key
!             integer, intent(in)              ::  num
!           end subroutine WriteKeyNum
!         \ev
!       \item
!         Replace tab-characters by spaces:
!         \bv
!           call goTab2Space( s )
!         \ev
!     \end{itemize}
!
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_string.F90} .
!
!
! !INTRODUCTION: GO_Date
!
!   General objects for date manipulation.
!
!   \subsection{Usage}
!
!     Access the entities of this module via the main module (prefered)
!     or directly:
!     \bv
!       use GO, only : ...
!       use GO_Date, only : ...
!     \ev
!
!
!   \subsection{Derived types}
!
!     A derived type is provided to store a date:
!     \bv
!       type(TDate)     ::  t
!     \ev
!     with fields:
!     \bv
!       character(len=4)   ::  calender     ! see 'calenders'
!       integer            ::  year, month, day, hour, min, sec, mili
!     \ev
!
!     A second type is defined for time steps:
!     \bv
!       type(TIncrDate)    ::  dt
!     \ev
!     with fields for days, hours, minutes, seconds, and mili-seconds.
!     The fields of an incremental data might have any value,
!     even negative. Note that for this type a value day=1
!     has the interpretation of 24 hours,
!     while it has the interpretation of 'first day', thus 0 hours,
!     for a regular date.
!
!
!   \subsection{Calenders}
!
!     A number of different calender types is supported:
!     \begin{itemize}
!       \item 'greg'
!         The Gregorian calender we are used to for a few ages now. \newline
!         Some years have a Februari 29. \newline
!         This is the default calender type; see also 'Defaults' below.
!       \item '366'
!         Every year has a Februari 29.
!       \item '365'
!         A year never has a  Februari 29.
!       \item '360'
!         Every month has 30 days. \newline
!         Use this calender if only operations on years and months are required.
!       \item 'wall'
!         Wall clock time, including time zone information.
!         Most operations are not implemented for this calender yet.
!     \end{itemize}
!
!
!
!   \subsection{Creating date structures}
!
!     To initialize a new date structure, a few routines are available.
!
!     Use routine \texttt{NewDate} to initialize some fields and to fill
!     the rest with zero's. If no calender is specified,
!     the default value \texttt{'greg'} is used (see section about defaults).
!     \bv
!       t = NewDate( year=2000, month=1, ..., calender='greg' )
!     \ev
!     Interface:
!     \bv
!       function NewDate( year, month, day, hour, min, sec, mili, zone, calender )
!         type(TDate)                              ::  NewDate
!         integer, intent(in), optional            ::  year, month, day
!         integer, intent(in), optional            ::  hour, min, sec, mili
!         integer, intent(in), optional            ::  zone
!         character(len=*), intent(in), optional   ::  calender
!       end function NewDate
!     \ev
!     Use the specific order of the optional fields to quickly set the
!     'largest' values only, for example the date:
!     \bv
!       t = NewDate( 2000, 1, 2 )
!     \ev
!     A special funtion creates a date that represents 'any' time;
!     useful to specify that some data is constant in time and thus
!     valid for any time. A special inquiry function is available to
!     test wether a date is any:
!     \bv
!       t = AnyDate()
!       l = IsAnyDate(t)
!     \ev
!
!     Use routine \texttt{IncrDate} to create an incremental date:
!     \bv
!       t = IncrDate( day=2, hour=12 )
!     \ev
!     Interface:
!     \bv
!       function IncrDate( day, hour, min, sec, mili )
!         type(TIncrDate)                  ::  IncrDate
!         integer, intent(in), optional    ::  day, hour, min, sec, mili
!       end function IncrDate
!     \ev
!
!     Fill the time from the system clock in a date structure;
!     the result is of calender type \texttt{'wall'}:
!     \bv
!       t = SystemDate()
!     \ev
!
!
!   \subsection{Changing date fields}
!
!     Use the \texttt{Set} routine to set some specific fields of
!     a date structure:
!     \bv
!       call Set( t [,year=2000] [,month=1] [,day=2]          &
!                   [,hour=0] [,min=0] [,sec=0] [,mili=0]     &
!                   [,zone=0] [calender='greg'] )
!     \ev
!     Use the specific order of the optional fields to quickly set the
!     'largest' values only, for example the date:
!     \bv
!       call Set( t, 2000, 1, 2 )
!     \ev
!
!     Normalize hours to 0,..,23, minutes to 0,..,59, etc:
!     \bv
!      call Normalize( t )
!     \ev
!
!     To check if all fields are consistent with eachother, use:
!     \bv
!       call Check( t )
!     \ev
!
!     Similar routines are implemented for date increments.
!
!
!   \subsection{Extraction of date fields}
!
!     Use the \texttt{Get} routine to extract some specific fields of
!     a date structure:
!     \bv
!       call Get( t [,year=year] [,month=month] [,day=day]              &
!                   [,hour=hour] [,min=min] [,sec=sec] [,mili=mili]     &
!                   [,zone=zone] [calender=calender] )
!     \ev
!     Use the specific order of the optional fields to quickly extract the
!     'largest' values only, for example the date:
!     \bv
!       call Get( t, year, month, day )
!     \ev
!
!     A similar routine is implemented for date increments.
!
!
!   \subsection{Inquiry functions}
!
!     A few inquiry functions are provided.
!
!     The logical function \texttt{LeapYear} tells you if the year
!     has a Februari 29 :
!     \bv
!       l = LeapYear( t )
!     \ev
!
!     Two integer functions are provided to count the total number
!     of days in a month or a year:
!     \bv
!       i = Days_in_Month( t )
!       i = Days_in_Year( t )
!     \ev
!
!     An integer function is provided to return the day number,
!     counting from 1 (Januari 1) to 360, 365, or 366 (last of December):
!     \bv
!       i = DayNumber( t )
!     \ev
!
!
!   \subsection{Operators}
!
!     Operator \texttt{+} is redefined to be able to add two
!     dates to each other.
!     Both should be of same calender type, unless one is an increment:
!     \bv
!        to =  t1 +  t2
!        to =  t  + dt
!       dto = dt1 + dt2
!     \ev
!
!     Operator \texttt{-} is redefined to substract a date from another.
!     Both should be of same calender type, unless the substracted date
!     (\texttt{t2}) is an increment:
!     \bv
!        to =  t1 -  t2
!        to =  t  - dt
!       dto = dt1 - dt2
!     \ev
!
!     Operator \texttt{*} has been redefined to multiply a date increment
!     with a real or integer number.
!     If necessary, a remaining fraction of miliseconds is rounded
!     to the nearest integer.
!     \bv
!       dto = dt * 2
!       dto = 2 * dt
!       dto = dt * 1.5
!       dto = 1.5 * dt
!     \ev
!
!     Operator \texttt{/} is redefined to devide a date incrment by a real
!     or an integer number.
!     If necessary, a remaining fraction of miliseconds is rounded
!     to the nearest integer.
!     \bv
!       dto = dt / 2
!       dto = dt / 1.5
!     \ev
!
!     Logical operators are defined to compare two dates with eachother;
!     both should be of same calender type:
!     \bv
!       t1 == t2
!       t1 >= t2
!       t1 >  t2
!       t1 <= t2
!       t1 <  t2
!     \ev
!
!
!   \subsection{Summation routines}
!
!     The total number in a certain unit is returned by \texttt{rTotal}
!     (real value) or \texttt{iTotal} (integer value, error if the sum could
!     only be expressed as fraction).
!     Currently supported units are \texttt{'year'}, \texttt{'month'}, \texttt{'day'},
!     \texttt{'hour'}, \texttt{'min'}, \texttt{'sec'}, and \texttt{'mili'}.
!     If the total number is not wel defined for a certain date
!     (how to assign a fraction of years to the date of today?),
!     an error message is produced.
!     Date increments are supported too.
!     \bv
!       r = rTotal( t, 'year'|'month'|... )
!       i = iTotal( t, 'year'|'month'|... )
!     \ev
!
!
!   \subsection{Temperal interpolation}
!
!     A linear interpolation in time is represented by:
!     \bv
!       f(t) = alfa1 * f(t1) + alfa2 * f(t2)
!     \ev
!     Given the dates t, t1, and t2, the fractions alfa1 and alfa2
!     are set by the following routine:
!     \bv
!       InterpolFractions( t, t1, t2, alfa1, alfa2, status )
!         type(TDate), intent(in)    ::  t
!         type(TDate), intent(in)    ::  t1
!         type(TDate), intent(in)    ::  t2
!         real, intent(out)          ::  alfa1
!         real, intent(out)          ::  alfa2
!         integer, intent(out)       ::  status
!     \ev
!
!
!   \subsection{Output}
!
!     To obtain a pretty formatted print of the value of a date,
!     the 'Pretty' routine is provided. Output differs based on
!     the calender type. Also implemented for date increments.
!     \bv
!       character(len=36) function Pretty( t )
!         type(TDate), intent(in)    ::  t
!
!       calender                      output
!       ----------------------------  ----------------------------------------
!       (date)
!         'wall'                       1:23:45:678 03 feb 2001 (GMT+02:00)
!         'greg', '366', '365', '360' 2001/02/03  1:23:45:678
!       (date increment)
!                                      2 days  1:23:45:678
!     \ev
!
!     Two routines are provided to write messages including a date to
!     the '\texttt{gol}' buffer from the '\texttt{GO_Print}' library.
!     Add a call to '\texttt{goPr}' or '\texttt{goErr}' to actually
!     display the message. Example:
!     \bv
!       call wrtgol( 'time  : ', t ); call goPr
!       call wrtgol( 'range : ', t1, ' - ', 't2' ); call goPr
!     \ev
!     provides:
!     \bv
!       [00] time  : 2001/01/01 00:00:00:000
!       [00] range : 2001/01/01 00:00:00:000 - 2001/01/01 03:00:00:000
!     \ev
!
!   \subsection{Defaults}
!
!     For setting some default values, the subroutine 'go_DateDefaults'
!     is available. All arguments are optional.
!     Yet, only the calender type might be set.
!     \bv
!       call goDateDefaults( [calender='greg'] )
!     \ev
!
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_date.F90} .
!
!
! !INTRODUCTION: GO_File
!
!   Selects free file units. \newline
!   Read lines from a commented text file, skipping comments and empty lines.
!
!   \subsection{Usage}
!
!     Access the entities of this module via the main module (prefered)
!     or directly:
!     \bv
!       use GO, only : ...
!       use GO_File, only : ...
!     \ev
!
!   \subsection{Free file unit}
!
!     Use the following routine to select a free (unopened) file unit:
!     \bv
!       call goGetFU( fu, status )
!         integer, intent(inout)    ::  fu
!         integer, intent(out)      ::  status
!     \ev
!     The routines searches for unopened file units with the range
!     specified by '\texttt{goFuRange}' from module '\texttt{GO_FU}',
!     excluding the standard units.
!     If all units are already opened, an error message is returned.
!
!   \subsection{Commented text files}
!
!     A commented text file might look like:
!     \bv
!       !
!       ! This is data file.
!       !
!
!       ! number of data values:
!       NDATA
!       3
!
!       ! data values:
!       DATA
!       1.0 2.0 3.0
!     \ev
!     To be able to read from this file without bothering about comments
!     and empty lines, a special file type is introduced:
!     \bv
!       type(TTextFile)    ::  file
!     \ev
!
!     To open the file, use:
!     \bv
!       subroutine Init( file, 'data.txt', iostat, status, comment )
!         type(TTextFile), intent(out)              ::  file
!         character(len=*), intent(in)              ::  filename
!         integer, intent(out)                      ::  iostat
!         character(len=*), intent(in), optional    ::  status
!         character(len=1), intent(in), optional    ::  comment
!     \ev
!     The 'iostat' argument is the 'usuall' status argument, thus on return
!     non-zero in case of errors.
!     The 'status' argument is the same as used by Fortran's '\texttt{open}' command
!     to specify wether a file should already exist (status='old'),
!     should not exist yet and will be created ('new') or might exist or not ('unknown').
!
!     The optional comment is a single character; lines in the file that
!     start with this character are skipped while reading.
!
!     To close the file, use:
!     \bv
!       subroutine file_Done( file, status )
!         type(TTextFile), intent(inout)    ::  file
!         integer, intent(out)              ::  status
!     \ev
!
!     To read one line, skipping empty and comment lines, use:
!     \bv
!       subroutine ReadLine( file, s, status  )
!         type(TTextFile), intent(inout)      ::  file
!         character(len=*), intent(out)       ::  s
!         integer, intent(out)                ::  status
!     \ev
!     Return status<0 means that no values are read since end of file is reached;
!     errors are indicated by status>0.
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_file.F90} .
!
!
!
! !INTRODUCTION: GO_Rc
!
!   Read settings from a resource file. \newline
!
!   \subsection{Resource files}
!
!     In the GO library, an 'rcfile' or 'resource' file is a text file
!     with settings for a program, with a format following the X-resource conventions:
!     \bv
!       ! This is an example resource file.
!       ! Use line of the form:
!       !
!       ! <name> : <value>
!       !
!       ! Some conventions:
!       !   * comment lines start with '!'
!       !   * logical values are 'T' or 'F'  (without quotes ...)
!       !   * character strings do not contain quotes
!
!       prog.n            :  20
!       input.file        :  /data/test.dat
!       prog.debug        :  T
!     \ev
!
!
!   \subsection{Derived types}
!
!     A type has been derived to get access to a rcfile:
!     \bv
!       type TrcFile
!     \ev
!
!
!   \subsection{Routines}
!
!     To open a rcfile, use:
!     \bv
!       subroutine Init( rcfile, fname, status )
!         type(TrcFile), intent(out)    ::  rcfile
!         character(len=*), intent(in)  ::  fname
!         integer, intent(out)          ::  status
!     \ev
!     To close it again, use:
!     \bv
!       subroutine Done( rcfile,  status )
!         type(TrcFile), intent(inout)    ::  rcfile
!         integer, intent(out)            ::  status
!     \ev
!
!     To read the value asigned to a key, use:
!     \bv
!       subroutine ReadRc( rcfile, key, x, status, default )
!         type(TrcFile), intent(in)                   ::  rcfile
!         character(len=*), intent(in)                ::  key
!         <type>, intent(out)                         ::  x
!         integer, intent(out)                        ::  status
!         <type>, intent(in), optional                ::  default
!     \ev
!     Current implemented types:
!     \bv
!       <type>  =  integer|real|character(len=*)|logical
!     \ev
!     If a key has not been found, an error messages is returned
!     unless a default is provided.
!
!   \subsection{See also}
!
!     The shell scripts '\texttt{go_readrc}' and '\texttt{go_pprc}' are available
!     to read from a rcfile withing shell scripts and to preprocess the rcfile
!     respectively.
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_rc.F90} .
!
!
!
! !INTRODUCTION: GO_System
!
!   \subsection{Description}
!
!     Interface to routines of some common system routines.
!     These routines are not part of the Fortran standard,
!     but almost always supplied.
!     There might be subtile compiler (vendor) specific differences however,
!     and therefor this module was made.
!     Preprocessor statements in select the compiler specific code.
!
!     Currently only implemented for:
!     \begin{itemize}
!       \item Intel Fortran Compiler (version 8)
!     \end{itemize}
!
!   \subsection{Usage}
!
!     Access the entities of this module via the main module (prefered)
!     or directly:
!     \bv
!       use GO, only : ...
!       use GO_System, only : ...
!     \ev
!
!
!   \subsection{Procedures}
!
!     \begin{itemize}
!       \item
!         Perform a system command and return exit status:
!         \bv
!           subroutine goSystem( command, status )
!             character(len=*), intent(in)     ::  command
!             integer, intent(inout)           ::  status
!           end subroutine goSystem
!         \ev
!       \item
!         Stop execution and set the return status:
!         \bv
!           subroutine goExit( status )
!             integer, intent(in)             ::  status
!           end subroutine goExit
!         \ev
!       \item
!         Return number of command line arguments:
!         \bv
!           integer function goArgC()
!           end function goArgC
!         \ev
!       \item
!         Return command line argument \texttt{nr} in character string \texttt{value}.
!         \bv
!           subroutine goGetArg( nr, value, status )
!             integer, intent(in)              ::  nr
!             character(len=*), intent(out)    ::  value
!             integer, intent(out)             ::  status
!           end subroutine goGetArg
!         \ev
!     \end{itemize}
!
!
!   \subsection{Hackers only}
!
!     The module is implemented in source file \texttt{go_system.F90} .
!
!
!EOI
!
!#######################################################################
!
#define TRACEBACK write (gol,'("ERROR in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
!#######################################################################


module GO

  use GO_FU

  use GO_Print

  use GO_System

  use GO_String

  use GO_Date
  use GO_Time

  use GO_File

  use GO_Rc

  use GO_Field

  !use goRestart

  !use goSave

  use GO_Timer, only : GO_Timer_Init, GO_Timer_Done
  use GO_Timer, only : GO_Timer_Def, GO_Timer_Start, GO_Timer_End

  use GO_Label, only : GO_Label_Init, GO_Label_Done, goLabel


  ! --- in/out ---------------------------------------

  public


  ! --- const ------------------------------------

  character(len=*), parameter, private   ::  mname = 'GO'


contains


  ! ====================================================================


  subroutine GO_Init( status )

    ! --- in/out -------------------------

    integer, intent(out)              ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Init'

    ! --- local --------------------------

    ! --- begin --------------------------

    ! ok
    status = 0

  end subroutine GO_Init


  ! ***


  subroutine GO_Done( status )

    ! --- in/out -------------------------

    integer, intent(out)       ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Done'

    ! --- local --------------------------

    ! --- begin --------------------------

    ! ok
    status = 0

  end subroutine GO_Done



end module GO


