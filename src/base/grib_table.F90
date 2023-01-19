!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!
!###############################################################################

module grib_table

  implicit none
  
  ! --- in/out -----------------------------------
  
  private
  
  public    ::  GetPid
  public    ::  GetPidName
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter  ::  mname = 'grib_table'
  


contains


  ! =============================================================


  integer function GetPid( table, key, status )
  
    use GO, only : gol, goErr
  
    ! --- in/out --------------------------------
  
    character(len=*), intent(in)        ::  table
    character(len=*), intent(in)        ::  key
    integer, intent(out)                ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GetPid'
    
    ! --- begin ----------------------------------
    
    select case ( table )
    
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ECMWF codes
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ec' )
      
        select case ( key )

          case ( 'CVL' , 'cvl'  ); GetPid = 027   ! low  vegetation cover
          case ( 'CVH' , 'cvh'  ); GetPid = 028   ! high vegetation cover
          case ( 'TVL' , 'tvl'  ); GetPid = 029   ! type of low  vegetation
          case ( 'TVH' , 'tvh'  ); GetPid = 030   ! type of high vegetation

          case ( 'CI'  , 'ci'   ); GetPid = 031   ! sea-ice cover

          case ( 'SWV1', 'swv1' ); GetPid = 039   ! volumetric soil water layer 1
          case ( 'SWV2', 'swv2' ); GetPid = 040   ! volumetric soil water layer 2
          case ( 'SWV3', 'swv3' ); GetPid = 041   ! volumetric soil water layer 3
          case ( 'SWV4', 'swv4' ); GetPid = 042   ! volumetric soil water layer 4

          case ( 'G10M', 'g10m' ); GetPid = 049   ! wind gust at 10m
          case ( '10FG', '10fg' ); GetPid = 049   ! wind gust at 10m

          case ( 'PV'           ); GetPid =  60   ! potential vorticity (K m2/kg/s)

          case ( 'MU'           ); GetPid = 104   !    updraught mass flux         kg m-2
          case ( 'MD'           ); GetPid = 105   !  downdraught mass flux         kg m-2
          case ( 'DU'           ); GetPid = 106   !    updraught detrainment rate  s m-1 
          case ( 'DD'           ); GetPid = 107   !  downdraught detrainment rate  s m-1 

          case ( 'K'            ); GetPid = 109   ! diffusion (m2/s)

          case ( 'U'            ); GetPid = 131   ! u wind
          case ( 'V'            ); GetPid = 132   ! v wind

          case ( 'Z'            ); GetPid = 129   ! geopotential, orography at surface
          case ( 'oro'          ); GetPid = 129   ! orography (geopotential at surface)
          case ( 'T'            ); GetPid = 130   ! temperature
          case ( 'Q'            ); GetPid = 133   ! specific humidity
          case ( 'SP'           ); GetPid = 134
          case ( 'W'            ); GetPid = 135
          case ( 'VO'           ); GetPid = 138

          case ( 'SLW' , 'slw'  ); GetPid = 140
          case ( 'SD'  , 'sd'   ); GetPid = 141
          case ( 'LSP' , 'lsp'  ); GetPid = 142
          case ( 'CP'  , 'cp'   ); GetPid = 143
          case ( 'SF'  , 'sf'   ); GetPid = 144
          case ( 'SSHF', 'sshf' ); GetPid = 146
          case ( 'SLHF', 'slhf' ); GetPid = 147  ! surface latent heat flux  (W m**-2 s)

          case ( 'LNSP', 'lnsp' ); GetPid = 152
          case ( 'D'   , 'd'    ); GetPid = 155
          case ( 'GH'  , 'gh'   ); GetPid = 156    ! geopotential height

          case ( 'BLH' , 'blh'  ); GetPid = 159    ! boundary layer height

          case ( 'U10M', 'u10m' ); GetPid = 165
          case ( 'V10M', 'v10m' ); GetPid = 166
          case ( 'T2M' , 't2m'  ); GetPid = 167
          case ( 'D2M' , 'd2m'  ); GetPid = 168

          case ( 'LSM' , 'lsm'  ); GetPid = 172   ! land/sea mask

          case ( 'z0m'          ); GetPid = 173
          case ( 'SR'  , 'sr'   ); GetPid = 173
          case ( 'AL'  , 'al'   ); GetPid = 174  ! Albedo   (0-1)

          case ( 'EWSS', 'ewss' ); GetPid = 180
          case ( 'NSSS', 'nsss' ); GetPid = 181

          case ( 'LSRH', 'lsrh' ); GetPid = 234  ! Log. Surf.Roughn. length for Heat

          case ( 'ustar'          ); GetPid = 250
          case ( 'Raero'          ); GetPid = 251

          case ( 'CLWC' , 'clwc'  ); GetPid = 246  ! cloud liquid water content
          case ( 'CIWC' , 'ciwc'  ); GetPid = 247  ! cloud ice water content
          case ( 'CC'   , 'cc'    ); GetPid = 248  ! cloud cover

          case ( 'SSR'  , 'ssr'   ); GetPid = 176  ! surface solar radiation
          case ( 'SRC'  , 'src'   ); GetPid = 198  ! skin reservoir content

          case default
            write (gol,'("unknown key `",a,"`")') trim(key); call goErr
            write (gol,'("  table : ",a)') trim(table); call goErr
            write (gol,'("in ",a)') rname; call goErr; status=1; return
        end select

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TM codes
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tm' )
      
        select case ( key )

          case ( 'PV'           ); GetPid =  60   ! potential vorticity (K m2/kg/s)

          case ( 'K'            ); GetPid = 109   ! diffusion (m2/s)

          case ( 'u'            ); GetPid = 131   ! u wind
          case ( 'v'            ); GetPid = 132   ! v wind

          case ( 'ps'           ); GetPid = 160   ! average surface pressure
          case ( 'pu'           ); GetPid = 158   ! u flux
          case ( 'pv'           ); GetPid = 159   ! v flux
          case ( 'pw'           ); GetPid = 157   ! w flux

          case ( 'Z'            ); GetPid = 129   ! geopotential, orography at surface

          case ( 'T'            ); GetPid = 130   ! temperature

          case ( 'Q'            ); GetPid = 133   ! specific humidity

          case ( 'GH'           ); GetPid = 156    ! geopotential height

          case ( 'eu'           ); GetPid = 161    ! entrainment updraft
          case ( 'du'           ); GetPid = 162    ! detrainment updraft
          case ( 'ed'           ); GetPid = 164    ! entrainment downdraft
          case ( 'dd'           ); GetPid = 165    ! detrainment downdraft

          case ( 'dk'           ); GetPid = 163    ! vertical diffusion coef.

          case ( 'clbas'        ); GetPid = 190    ! cloud base
          case ( 'cltop'        ); GetPid = 191    ! cloud top
          case ( 'cllfs'        ); GetPid = 192

          case ( 'Kz'           ); GetPid = 163

          case ( 'pblh'         ); GetPid = 159
          case ( 'pblh?'        ); GetPid = 203

          case ( 'lsp'          ); GetPid = 142 
          case ( 'cp'           ); GetPid = 143 
          case ( 'sf'           ); GetPid = 144 
          case ( 'sshf'         ); GetPid = 146 
          case ( 'slhf'         ); GetPid = 147  ! surface latent heat flux  (W m**-2 s)

          case ( 'slw'          ); GetPid = 140
          case ( 'sd'           ); GetPid = 141
          case ( 'T2M', 'T2m', 't2m' ); GetPid = 167
          case ( 'D2M', 'D2m', 'd2m' ); GetPid = 168
          case ( 'z0m', 'SR'    ); GetPid = 173
          case ( 'al'           ); GetPid = 174  ! Albedo    (0-1) 

          case ( 'ewss'         ); GetPid = 180                                             
          case ( 'nsss'         ); GetPid = 181                                             

          case ( 'ustar'        ); GetPid = 250                                             
          case ( 'Raero'        ); GetPid = 251                                             

          case ( 'CLWC', 'clwc' ); GetPid = 246  ! cloud liquid water content
          case ( 'CIWC', 'ciwc' ); GetPid = 247  ! cloud ice water content
          case ( 'CC'  , 'cc'   ); GetPid = 248  ! cloud cover
          case ( 'CCO' , 'cco'  ); GetPid = 249  ! overhead cloud cover
          case ( 'CCU' , 'ccu'  ); GetPid = 250  ! underfeet cloud cover

          case ( 'ssr'          ); GetPid = 176  ! surface solar radiation
          case ( 'src'          ); GetPid = 198  ! skin reservoir content

          case default
            write (gol,'("unknown key `",a,"`")') trim(key); call goErr
            write (gol,'("  table : ",a)') trim(table); call goErr
            write (gol,'("in ",a)') rname; call goErr; GetPid=-1; return
        end select

      case default
        write (gol,'("unknown table `",a,"`")') trim(table); call goErr
        write (gol,'("in ",a)') rname; call goErr; GetPid=-1; return
    end select

  end function GetPid

  
  ! =============================================================


  character(len=4) function GetPidName( table, pid )
  
    use GO, only : gol, goErr
  
    ! --- in/out ------------------------------
    
    character(len=*), intent(in)        ::  table
    integer, intent(in)                 ::  pid
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GetPidName'
    
    ! --- local ------------------------------
    
    character(len=4)      ::  res
    
    ! --- begin -------------------------------
    
    select case ( table )
    
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ECMWF codes
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'ec' )
      
        select case ( pid )
          case ( 109 )  ;  res = 'K   '     ! diffusion (m2/s); experimental
          case ( 129 )  ;  res = 'Z   '
          case ( 130 )  ;  res = 'T   '
          case ( 131 )  ;  res = 'U   '
          case ( 132 )  ;  res = 'V   '
          case ( 133 )  ;  res = 'Q   '
          case ( 134 )  ;  res = 'SP  '
          case ( 135 )  ;  res = 'W   '
          case ( 138 )  ;  res = 'VO  '
          case ( 141 )  ;  res = 'SD  '
          case ( 142 )  ;  res = 'LSP '
          case ( 143 )  ;  res = 'CP  '
          case ( 144 )  ;  res = 'SF  '
          case ( 146 )  ;  res = 'SSHF'
          case ( 147 )  ;  res = 'SLHF'
          case ( 152 )  ;  res = 'LNSP'
          case ( 155 )  ;  res = 'D   '
          case ( 156 )  ;  res = 'GH  '
          case ( 159 )  ;  res = 'BLH '
          case ( 165 )  ;  res = 'U10M'
          case ( 166 )  ;  res = 'V10M'
          case ( 167 )  ;  res = 'T2M '
          case ( 168 )  ;  res = 'D2M '
          case ( 172 )  ;  res = 'LSM '
          case ( 173 )  ;  res = 'SR  '
          case ( 174 )  ;  res = 'AL  '
          case ( 176 )  ;  res = 'SSR '
          case ( 180 )  ;  res = 'EWSS'
          case ( 181 )  ;  res = 'NSSS'
          case ( 198 )  ;  res = 'SRC '
          case ( 234 )  ;  res = 'LSRH'
          case ( 246 )  ;  res = 'CLWC'
          case ( 247 )  ;  res = 'CIWC'
          case ( 248 )  ;  res = 'CC  '
          case default
            write (res,'("p",i3.3)') pid
        end select
    

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TM codes
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      case ( 'tm' )
      
        select case ( pid )
          case ( 104 )  ;  res = 'MU  '     !    updraught mass flux         kg m-2
          case ( 105 )  ;  res = 'MD  '     !  downdraught mass flux         kg m-2
          case ( 106 )  ;  res = 'DU  '     !    updraught detrainment rate  s m-1 
          case ( 107 )  ;  res = 'DD  '     !  downdraught detrainment rate  s m-1 
          case ( 109 )  ;  res = 'K   '     ! diffusion (m2/s); experimental
          case ( 129 )  ;  res = 'Z   '
          case ( 130 )  ;  res = 'T   '
          case ( 131 )  ;  res = 'U   '
          case ( 132 )  ;  res = 'V   '
          case ( 133 )  ;  res = 'Q   '
          case ( 134 )  ;  res = 'SP  '
          case ( 135 )  ;  res = 'W   '
          case ( 138 )  ;  res = 'VO  '
          case ( 141 )  ;  res = 'SD  '
          case ( 142 )  ;  res = 'LSP '
          case ( 143 )  ;  res = 'CP  '
          case ( 144 )  ;  res = 'SF  '
          case ( 146 )  ;  res = 'SSHF'
          case ( 147 )  ;  res = 'SLHF'
          case ( 152 )  ;  res = 'LNSP'
          case ( 155 )  ;  res = 'D   '
          case ( 156 )  ;  res = 'zg  '      ! ECMWF : GH
          case ( 157 )  ;  res = 'pw  '
          case ( 158 )  ;  res = 'pu  '
          case ( 159 )  ;  res = 'pv  '
          case ( 160 )  ;  res = 'sp  '
          case ( 161 )  ;  res = 'eu  '
          case ( 162 )  ;  res = 'du  '
          case ( 163 )  ;  res = 'dk  '
          case ( 164 )  ;  res = 'ed  '
          case ( 165 )  ;  res = 'dd  '
          case ( 167 )  ;  res = 'T2M '
          case ( 168 )  ;  res = 'D2M '
          case ( 172 )  ;  res = 'LSM '
          case ( 173 )  ;  res = 'SR  '
          case ( 174 )  ;  res = 'AL  '
          case ( 176 )  ;  res = 'SSR '
          case ( 180 )  ;  res = 'EWSS'
          case ( 181 )  ;  res = 'NSSS'
          case ( 190 )  ;  res = 'clb '
          case ( 191 )  ;  res = 'clt '
          case ( 192 )  ;  res = 'clfs'
          case ( 198 )  ;  res = 'SRC '
          case ( 203 )  ;  res = 'blh '     ! ECMWF's BLH has code 159, but this is pv in TM ...
          case ( 234 )  ;  res = 'LSRH'
          case ( 246 )  ;  res = 'CLWC'
          case ( 247 )  ;  res = 'CIWC'
          case ( 248 )  ;  res = 'CC  '
          case ( 249 )  ;  res = 'cco '
          case ( 250 )  ;  res = 'ccu '
          case default
            write (res,'("p",i3.3)') pid
        end select
    
      case default
        write (gol,'("unknown table `",a,"`")') trim(table); call goErr
        write (gol,'("in ",a)') rname; call goErr; GetPidName='ERROR'; return
    end select

    GetPidName = res
    
  end function GetPidName
  

end module grib_table
