!======================================================================
!
!   Name:               utilities.f90
!
!======================================================================
module mo_netcdf
  use netcdf 
  character(len=5) :: ipath = 'input'    ! input directory
  character(len=3+3) :: foj = 'foj.nc'        ! obs jacobian file
  character(len=3+3) :: ftj = 'ftj.nc'        ! target jacobian file
  character(len=2+3) :: fe = 'fe.nc'          ! prior emissions file
  character(len=7+3) :: fesigma = 'fesigma.nc'          ! prior emissions file
  character(len=2+3) :: fc = 'fc.nc'          ! simulated concentration at sites
  character(len=6+3) :: fcpost = 'fcpost.nc'          ! simulated concentrations prior+post
  character(len=6+3) :: fepost = 'fepost.nc'          ! posterior emissions file
  integer ( kind = 4 ) :: nobsday, nsta, ng, ntgt
  real ( kind = 8 ) :: kg2mt=1.e-9 ! unit conversion factor
end module mo_netcdf

subroutine rc_handler(status,erc)
! - check NetCDF error status,
!   optionally stop running program
  use netcdf
  implicit none
  integer, intent(in) :: status
  logical, intent(in) :: erc
  if(status /= nf90_noerr) then 
     write(*,'(a,a)') 'ERROR::',trim(nf90_strerror(status))
     if( erc ) then
        stop 1
     end if
  end if
end subroutine rc_handler

subroutine wfmat(unit,fname,dim1,dim2,mat,description)
  implicit none
  integer ( kind = 4 ), intent(in) :: unit, dim1, dim2
  real ( kind = 8 ), intent(in) :: mat(dim1,dim2)
  character(len=*), intent(in) :: description
  character(len=*), intent(in) :: fname
  
! local
  integer ( kind = 4 ) :: i,j

  if (unit.eq.6) then
     write ( unit, '(a)' ) ' '
     write ( unit, '(a)' ) description
     write ( unit, '(a)' ) 'output limited to 5 rows and 5 columns'
  else
     open (unit=unit,file=fname, form='formatted', action='write')
  endif
  write ( unit, '(a9,1x,5(a19,i1))') 'row_#', ('               col_',i,i=1,min(5,dim2))
  do j = 1, min(5,dim1)
     write ( unit, '(1i9,1x,5f20.8)') j, (mat(j,i),i=1,min(5,dim2))
  enddo
  if (unit.eq.6) then
     write ( unit, '(a10,1x,f20.8,a12,2i10)') 'maximum = ', maxval(mat) , ' at indices ', maxloc(mat)
     write ( unit, '(a10,1x,f20.8,a12,2i10)') 'minimum = ', minval(mat) , ' at indices ', minloc(mat)
     write ( unit, '(a)' ) ' '
  else
     close (unit=unit)
  endif
end subroutine wfmat
subroutine setdim(n,m,k)
  ! specifies dimension of the problem
  implicit none
  integer ( kind = 4 ), intent(out)    :: n     ! dim of control space
  integer ( kind = 4 ), intent(out)    :: m     ! dim of obs space
  integer ( kind = 4 ), intent(out)    :: k     ! dim of target space

  n = 60*45*12 ! global control
!  n = 20*45*12 ! a small value for testing
!  n = 1440+432+288 ! monthly for 1 month
  n = n + 30*26*52 ! European domain
  n = n + 18*16*365 ! Zoom domain
  m = 20*12 ! global obs
!  m = m + 30*10 ! European domain obs
!  m = m*2
!  m = m + 8*365 ! Zoom domain obs
!  m = 31*4 ! for 4 stations in zoom domain
!  m = m + 1*5 ! for 5 stations in the rest of the world
!  m = 20 ! a small value for testing
!  m = m + 30*365 ! Zoom domain obs
  k = 3 ! could start with emissions total (over sectors and target period)
        ! for NL, D, and CH, can later be changed to sectoral
  k = k+1 ! and the global total

  ! I think for the first inversion, we could use the Feb Jacobian for Cabauw,
  ! i.e. 31 observations and
  ! monthly emissions outside and inside target region, i.e. 2160 (1440+432+288) and
  ! targets as suggested above
  ! in the next step we can then move to daily emissions inside target region
  
end subroutine setdim

subroutine setdim2(n,m,k)
  ! specifies dimension of the problem
  use mo_netcdf
  implicit none
  ! args
  integer ( kind = 4 ), intent(out)    :: n     ! dim of control space
  integer ( kind = 4 ), intent(out)    :: m     ! dim of obs space
  integer ( kind = 4 ), intent(out)    :: k     ! dim of target space
  ! local
  integer ( kind = 4 ) :: rc, ncid, dimid

  ! first obs jac, get ncid
  rc = nf90_open(foj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  ! get dimid ng
  rc = nf90_inq_dimid(ncid, "ng", dimid); call rc_handler(rc, .true.)
  ! get ng
  rc = nf90_inquire_dimension(ncid, dimid, len=ng)
  print*, 'ng = ', ng
  ! get dimid nsta
  rc = nf90_inq_dimid(ncid, "nsta", dimid); call rc_handler(rc, .true.)
  ! get nsta
  rc = nf90_inquire_dimension(ncid, dimid, len=nsta)
  print*, 'nsta = ', nsta
  ! get dimid nobsday
  rc = nf90_inq_dimid(ncid, "nobsday", dimid); call rc_handler(rc, .true.)
  ! get nobsday
  rc = nf90_inquire_dimension(ncid, dimid, len=nobsday)
  rc = nf90_close(ncid)
  print*, 'nobsday = ', nobsday
  ! now target jac, get ncid
  rc = nf90_open(ftj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  ! get dimid ntgt
  rc = nf90_inq_dimid(ncid, "ntgt", dimid); call rc_handler(rc, .true.)
  ! get ntgt
  rc = nf90_inquire_dimension(ncid, dimid, len=ntgt)
  rc = nf90_close(ncid)
  print*, 'ntgt = ', ntgt
  ! set outputs
  n=ng
  m=nsta*nobsday
  k=ntgt
end subroutine setdim2

subroutine setval (n, m, k, x0, sx0, xp, y, sy, yinit, oj, tj, twin)
  ! specifies values of required inputs
  implicit none
  integer ( kind = 4 ), intent(in)    :: n     ! dim of control space
  integer ( kind = 4 ), intent(in)    :: m     ! dim of obs space
  integer ( kind = 4 ), intent(in)    :: k     ! dim of target space
  real ( kind = 8 ), intent(out) :: x0(n), sx0(n) ! prior and uncertainty in kgCH4/m^2/s
  real ( kind = 8 ), intent(out) :: xp(n) ! perturbed prior for identical twin
  real ( kind = 8 ), intent(out) :: y(m), sy(m) ! obs contribution from emissions
                                                ! and uncertainty in ppb
  real ( kind = 8 ), intent(out) :: yinit(m) ! simulated from initial [ppb]
  real ( kind = 8 ), intent(out) :: oj(m,n) ! obs Jacobian in units of sy/sx0
  real ( kind = 8 ), intent(out) :: tj(k,n) ! target Jacobian in units kgCH4/sx0
  logical, intent(out) :: twin ! setup for twin experiment
! local
  integer ( kind = 4 ) :: i,j,l
  real ( kind = 8 ) :: pert=1.1_8 ! perturbation for identical twin

  x0 = 2._8; sx0 = 1_8 ! dummy settings
  y = 4_8; sy = 10._8   ! dummy settings, from y subtract simulated response to initial concentration
  ! dummy setting
  do i = 1, n
     do j = 1, m
        oj(j,i) = i+j
     enddo
  enddo
  ! set targets
  tj=0.
  do l = 1, k
     tj(l,1) = 1._8
  enddo
  tj(1,2) = 1._8
  tj(1,3) = 1._8
  tj(3,3) = 1._8
  tj(k,:) = 1._8
  ! so tj =
  !
  ! (1,1,1,0, ... , 0)
  ! (1,0,0,0, ... , 0)
  ! (1,0,1,0, ... , 0)
  ! ...
  ! (1,1,1,1, ... , 1)
  ! 
  ! normalisation obs jac
  do i = 1, n
     do j = 1, m
        oj(j,i) = oj(j,i)*sx0(i)/sy(j)
     enddo
  enddo  
  ! normalisation target jac
  do i = 1, n
     do l = 1, k
        tj(l,i) = tj(l,i)*sx0(i)
     enddo
  enddo
  ! normalisation prior and obs
  x0 = x0/sx0
  y = y/sy
  ! identical twin experiment with known truth
  if (twin) then
     xp = x0*pert
     y = matmul (oj,xp)
  else
     xp = 0.
  endif
end subroutine setval

subroutine setval2 (n, m, k, x0, sx0, xp, y, sy, yinit, oj, tj, twin)
  ! specifies values of required inputs
  use mo_netcdf
  implicit none
  integer ( kind = 4 ), intent(in)    :: n     ! dim of control space
  integer ( kind = 4 ), intent(in)    :: m     ! dim of obs space
  integer ( kind = 4 ), intent(in)    :: k     ! dim of target space
  real ( kind = 8 ), intent(out) :: x0(n), sx0(n) ! prior and uncertainty in MtCH4/cell
  real ( kind = 8 ), intent(out) :: xp(n) ! perturbed prior for identical twin
  real ( kind = 8 ), intent(out) :: y(m), sy(m) ! obs contribution from emissions
                                                ! and uncertainty in ppb
  real ( kind = 8 ), intent(out) :: yinit(m) ! simulated from initial [ppb]
  real ( kind = 8 ), intent(out) :: oj(m,n) ! obs Jacobian in units of sy/sx0
  real ( kind = 8 ), intent(out) :: tj(k,n) ! target Jacobian in units MtCH4/sx0
  logical, intent(out) :: twin ! setup for twin experiment
! local
  integer ( kind = 4 ) :: i,j,l,ista
  real ( kind = 8 ) :: pert=1.1_8 ! perturbation for identical twin
  real ( kind = 8 ) :: minsx0 ! floor for prior uncertainty
  integer ( kind = 4 ) :: rc, ncid, varid, no_fill
  real ( kind = 8 ) :: fill_value ! 
  real ( kind = 8 ) :: ojin(ng, nsta, nobsday) ! input obs Jacobian [ppb/(kgCH4/cell)]
  real ( kind = 8 ) :: yin(nsta, nobsday) ! input obs [ppb]
  real ( kind = 8 ) :: tjin(ng, ntgt) ! input target Jacobian [kgCH4/(kgCH4/cell)]
  logical :: mask(nsta,nobsday) ! to be set false for fill values
  logical :: complete(nsta) ! to be set false if any fill value
  
  x0 = 2._8; sx0 = 1_8 ! dummy settings
  y = 4_8; sy = 10._8   ! dummy settings, from y subtract simulated response to initial concentration
  ! first emissions
  rc = nf90_open(fe,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'emission', varid); call rc_handler(rc, .true.)
  rc = nf90_get_var(ncid, varid, x0); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  call wfmat(6,'dummy.dat',1,n,x0,'x0') ! this is just for development
  x0 = x0*kg2mt
  call wfmat(6,'dummy.dat',1,n,x0,'x0 in MtCH4') ! this is just for development
  ! ad hoc setting for emission uncertainty
  rc = nf90_open(fesigma,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'emission', varid); call rc_handler(rc, .true.)
  rc = nf90_get_var(ncid, varid, sx0); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  call wfmat(6,'dummy.dat',1,n,sx0,'sx0') ! this is just for development
  sx0 = sx0*kg2mt
  call wfmat(6,'dummy.dat',1,n,x0,'sx0 in MtCH4') ! this is just for development
  print*, 'sum(sx0), 0.01*sum(sx0)/n = ', sum(sx0), 0.01*sum(sx0)/n
  minsx0=0.01*sum(sx0)/n
  where (abs(sx0).gt.minsx0)
     sx0=3*abs(sx0)
  elsewhere
     sx0=minsx0
  endwhere
  print*, 'minval(sx0) = ', minval(sx0)
  ! now concentration, observed and simulated from initial
  rc = nf90_open(foj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'obs', varid); call rc_handler(rc, .true.)
  !  print*, 'ncid, varid = ', ncid, varid
  rc = nf90_get_var(ncid, varid, yin); call rc_handler(rc, .true.)
! defining mask to indicated missing values
  mask=.true.
  where (yin.ne.yin)
     mask=.false.
  endwhere
! count sites without missing values
  write(11,'(a10,1x,a10)') 'site#', 'complete'
  do ista=1,nsta
     write(11,'(i10,1x,31l10)') ista,all(mask(ista,:))
     complete(ista) = all(mask(ista,:))
  enddo
  write(11,'(a,i10)') '# of sites with complete obs recored = ', count(complete)
  call flush(11)
  print*, '# of valid obs = ', count(mask), ' from ', nobsday*nsta
  write(10,'(a10,1x,31i10)') 'site#', (i,i=1,31)
  do ista=1,nsta
     write(10,'(i10,1x,31f10.1)') ista,yin(ista,:)
     write(10,'(i10,1x,31l10)') ista,mask(ista,:)
  enddo
  ! FILL_VALUE=99.
  ! rc = NF90_INQ_VAR_FILL(ncid, varid, no_fill, FILL_VALUE); call rc_handler(rc, .true.)
  ! print*, ' no_fill, FILL_VALUE =', no_fill, FILL_VALUE
  ! call wts('output/ytest1.dat',nobsday,yin(1,:))
  ! if (nsta.gt.1) call wts('output/ytest2.dat',nobsday,yin(2,:))
  ! replace missing values by -1 
  y = reshape(yin,(/nobsday*nsta/))
!  stop
  call wfmat(6,'dummy.dat',1,m,y,'y') ! this is just for development
  rc = nf90_inq_varid(ncid, 'iniconc', varid); call rc_handler(rc, .true.)
!  print*, 'ncid, varid = ', ncid, varid
  rc = nf90_get_var(ncid, varid, yin); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  yinit = reshape(yin,(/nobsday*nsta/))
  call wfmat(6,'dummy.dat',1,m,yinit,'yinit') ! this is just for development
  y=y-yinit 
  call wfmat(6,'dummy.dat',1,m,y,'y-yinit') ! this is just for development
  ! now obs jac
  ! dummy setting
  do i = 1, n
     do j = 1, m
        oj(j,i) = i+j
     enddo
  enddo
  ! now reading
  rc = nf90_open(foj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'obs_jacobian', varid); call rc_handler(rc, .true.)
!  print*, 'ncid, varid = ', ncid, varid
  rc = nf90_get_var(ncid, varid, ojin); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  oj = transpose(reshape(ojin,(/ng,nobsday*nsta/)))
  print*, 'shape oj = ', shape(oj)
  call wfmat(6,'dummy.dat',m,n,oj,'oj') ! this is just for development
  oj=oj/kg2mt
  call wfmat(6,'dummy.dat',m,n,oj,'oj in ppp/MtCH4/cell') ! this is just for development
  ! set targets
  tj=0.
  do l = 1, k
     tj(l,1) = 1._8
  enddo
  tj(1,2) = 1._8
  tj(1,3) = 1._8
  tj(k,:) = 1._8
  ! so tj =
  !
  ! (1,1,1,0, ... , 0)
  ! (1,0,0,0, ... , 0)
  ! ...
  ! (1,1,1,1, ... , 1)
  ! 
  ! now reading
  rc = nf90_open(ftj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'tgt_jacobian', varid); call rc_handler(rc, .true.)
!  print*, 'ncid, varid = ', ncid, varid
  rc = nf90_get_var(ncid, varid, tjin); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  tj = transpose(tjin)
  print*, 'shape tj = ', shape(tj)
  call wfmat(6,'dummy.dat',k,n,tj,'tj') ! this is just for development
  ! normalisation obs jac
  do i = 1, n
     do j = 1, m
        oj(j,i) = oj(j,i)*sx0(i)/sy(j)
     enddo
  enddo  
  call wfmat(6,'dummy.dat',m,n,oj,'oj scaled') ! this is just for development
  ! normalisation target jac
  do i = 1, n
     do l = 1, k
        tj(l,i) = tj(l,i)*sx0(i)
     enddo
  enddo
  call wfmat(6,'dummy.dat',k,n,tj,'tj scaled') ! this is just for development
  ! normalisation prior and obs
  x0 = x0/sx0
  call wfmat(6,'dummy.dat',1,n,x0,'x0 scaled') ! this is just for development
  y = y/sy
  ! identical twin experiment with known truth
  if (twin) then
     xp = x0*pert
     call wfmat(6,'dummy.dat',1,n,xp,'xp') ! this is just for development
     y = matmul (oj,xp)
     call wfmat(6,'dummy.dat',1,m,y,'y twin') ! this is just for development
  else
     xp = 0.
  endif
end subroutine setval2

subroutine setvalf (n, m, k, x0, yinit, oj)
  ! specifies values of required inputs
  use mo_netcdf
  implicit none
  integer ( kind = 4 ), intent(in)    :: n     ! dim of control space
  integer ( kind = 4 ), intent(in)    :: m     ! dim of obs space
  integer ( kind = 4 ), intent(in)    :: k     ! dim of target space
  real ( kind = 8 ), intent(out) :: x0(n) ! emission field [kgCH4/cell]
  real ( kind = 8 ), intent(out) :: yinit(m) ! simulated from initial [ppb]
  real ( kind = 8 ), intent(out) :: oj(m,n) ! obs Jacobian 
! local
  integer ( kind = 4 ) :: rc, ncid, varid
  real ( kind = 8 ) :: ojin(ng, nsta, nobsday) ! input obs Jacobian [ppb/(kgCH4/cell)]
  real ( kind = 8 ) :: yin(nsta, nobsday) ! input obs [ppb]

  ! first emissions
  rc = nf90_open(fe,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'emission', varid); call rc_handler(rc, .true.)
  rc = nf90_get_var(ncid, varid, x0); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  call wfmat(6,'dummy.dat',1,n,x0,'x0') ! this is just for development
  ! now concentration, simulated from initial
  rc = nf90_open(foj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'iniconc', varid); call rc_handler(rc, .true.)
  rc = nf90_get_var(ncid, varid, yin); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  yinit = reshape(yin,(/nobsday*nsta/))
  call wfmat(6,'dummy.dat',1,m,yinit,'yinit') ! this is just for development
  ! now obs jac
  rc = nf90_open(foj,nf90_nowrite, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'obs_jacobian', varid); call rc_handler(rc, .true.)
  rc = nf90_get_var(ncid, varid, ojin); call rc_handler(rc, .true.)
  rc = nf90_close(ncid)
  oj = transpose(reshape(ojin,(/ng,nobsday*nsta/)))
end subroutine setvalf

subroutine wemission(n,x)
  ! write emission in NetCDF 
  use mo_netcdf
  implicit none
  integer ( kind = 4 ), intent(in)    :: n     ! dim of control space
  real ( kind = 8 ), intent(in) :: x(n) ! posterior emission in MtCH4/cell
  ! local
  integer ( kind = 4 ) :: rc, ncid, varid
  ! open file
  rc = nf90_open(fepost,nf90_write, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'emission', varid); call rc_handler(rc, .true.)
  rc = nf90_put_var(ncid, varid, x/kg2mt) ! back to kgCH4/cell
  rc = nf90_close(ncid)
end subroutine wemission

subroutine wtarget(k,t0, deltat,ct0, ct)
  ! specifies posterior targets and their uncertainty
  implicit none
  integer ( kind = 4 ), intent(in)    :: k     ! dim of target space
  real ( kind = 8 ), intent(in)       :: t0(k) ! prior targets
  real ( kind = 8 ), intent(in)       :: deltat(k) ! delta from inversion
  real ( kind = 8 ), intent(in)       :: ct0(k,k)  ! covar of prior target uncertainty
  real ( kind = 8 ), intent(in)       :: ct(k,k)   ! covar of posterior target uncertainty
  ! local
  integer ( kind = 4 )    :: l    ! counter

  ! output targets and their uncertainty such that they can be displayed by GUI
  ! and later also ensure output in common data format
  ! (agreed between AVENGERS, PARIS, and EYE-CLIMA)
  call wfmat(6,'dummy.dat',k,1,t0,'prior target') ! this is just for development
  call wfmat(6,'dummy.dat',k,1,t0+deltat,'posterior target') ! this is just for development
  call wfmat(6,'dummy.dat',k,k,ct0,'prior target uncertainty') ! this is just for development
  call wfmat(6,'dummy.dat',k,k,ct,'posterior target uncertainty') ! this is just for development
  call wfmat(10,'output/t0.dat',k,1,t0,'prior target') ! this is just for development
  call wfmat(10,'output/t.dat',k,1,t0+deltat,'posterior target') ! this is just for development
  call wfmat(10,'output/ct0.dat',k,k,ct0,'prior target uncertainty') ! this is just for development
  call wfmat(10,'output/ct.dat',k,k,ct,'posterior target uncertainty') ! this is just for development
  write (6, '(1a9,1x,4a12)') 'target #', 'prior', 'posterior', 'prior sigma', 'post. sigma'
  do l=1,k
     write (6, '(1i9,1x,4f12.6)') l , t0(l), t0(l)+deltat(l), sqrt(ct0(l,l)), sqrt(ct(l,l))
  enddo

end subroutine wtarget

subroutine wconc(m,y)
  ! write concentration in NetCDF 
  use mo_netcdf
  implicit none
  integer ( kind = 4 ), intent(in)    :: m   ! dim of obs
  real ( kind = 8 ), intent(in) :: y(m) ! simulated concentration in ppb
  ! local
  integer ( kind = 4 ) :: rc, ncid, varid
  real ( kind = 8 ) :: yout(nsta,nobsday) ! for writing onto netcdf
  ! open file
  rc = nf90_open(fc,nf90_write, ncid); call rc_handler(rc, .true.)
  rc = nf90_inq_varid(ncid, 'obs', varid); call rc_handler(rc, .true.)
  yout = reshape(y,(/nsta,nobsday/))
  rc = nf90_put_var(ncid, varid, yout)
  rc = nf90_rename_var(ncid, varid, 'conc')
  rc = nf90_close(ncid)
end subroutine wconc

subroutine wconcpost(m,y0,y,yobs)
  ! write concentration in NetCDF 
  use mo_netcdf
  implicit none
  integer ( kind = 4 ), intent(in)    :: m   ! dim of obs
  real ( kind = 8 ), intent(in) :: y0(m) ! concentration simulated from prior [ppb]
  real ( kind = 8 ), intent(in) :: y(m) ! concentration simulated from posterior [ppb]
  real ( kind = 8 ), intent(in) :: yobs(m) ! obs [ppb]
  ! local
  integer ( kind = 4 ) :: rc, ncid, varid, nobsdaydimid, nstadimid, yvarid, y0varid
  real ( kind = 8 ) :: yout(nsta,nobsday) ! for writing onto netcdf
  print*, 'in wconc: '
  print*, 'RMSE cprior = ', sqrt(sum(((y0-yobs))**2)/m)
  print*, 'RMSE cpost = ', sqrt(sum(((y-yobs))**2)/m)
  ! open file
  ! print*, 'wconcpost@'//trim(fcpost)
  rc = nf90_open(fcpost,nf90_write, ncid); call rc_handler(rc, .true.)
  ! get dimids
  rc = nf90_inq_dimid(ncid, "nobsday", nobsdaydimid); call rc_handler(rc, .true.)
  rc = nf90_inq_dimid(ncid, "nsta", nstadimid); call rc_handler(rc, .true.)
  rc = nf90_enddef(ncid)
  ! define variables
  !-- ATTENTION: station is the faster iterating dimension and must
  !              become first dimension in *Fortran*
  rc = nf90_def_var(ncid, "cpost", nf90_real8, (/ nstadimid, nobsdaydimid /), yvarid)
  call rc_handler(rc, .true.)
  rc = nf90_def_var(ncid, "cprior", nf90_real8, (/ nstadimid, nobsdaydimid /), y0varid)
  call rc_handler(rc, .true.)
  print*, 'conc defined'; call flush(6)
  ! set attributes
  rc = nf90_put_att(ncid, yvarid, "long_name", "concentration simulated from posterior")
  rc = nf90_put_att(ncid, yvarid, "units", "ppb")
  rc = nf90_put_att(ncid, y0varid, "long_name", "concentration simulated from prior")
  rc = nf90_put_att(ncid, y0varid, "units", "ppb")
  ! write variables
  yout = reshape(y,(/nsta,nobsday/))
  ! print*, 'wconcpost: writing posterior to yvarid=',yvarid, 'shape(transpose(yout))=', &
  !      shape(transpose(yout))
  rc = nf90_put_var(ncid, yvarid, yout)
  call rc_handler(rc, .true.) !-- check writing data went properly 
  yout = reshape(y0,(/nsta,nobsday/))
  rc = nf90_put_var(ncid, y0varid, yout)
  call rc_handler(rc, .true.) !-- check writing data went properly 
  ! close file
  rc = nf90_close(ncid)
end subroutine wconcpost

! status = nf90_create("foo.nc", nf90_noclobber, ncid)
! if(status /= nf90_noerr) call handle_error(status)
! ...
! ! Define the dimensions
! status = nf90_def_dim(ncid, "lat", 5, latdimid)
! if(status /= nf90_noerr) call handle_error(status)
! status = nf90_def_dim(ncid, "lon", 10, londimid)
! if(status /= nf90_noerr) call handle_error(status)
! status = nf90_def_dim(ncid, "time", nf90_unlimited, timedimid)
! if(status /= nf90_noerr) call handle_error(status)
! ...
! ! Define the variable
! status = nf90_def_var(ncid, "rh", nf90_double, &
!                       (/ londimid, latdimid, timedimid /), rhvarid)
! if(status /= nf90_noerr) call handle_error(status)

subroutine wts(fname, m, y)
! writing observational data stream
  implicit none
  integer, intent(in) :: m
  real(kind=8), intent(in) :: y(m)
  character(*) :: fname
  integer :: i
  logical :: debug=.false.

  if (debug) print*, 'Opening file ',fname,' now'
  open(1,file=fname,action="write",form='formatted')
  write(1,'(a,a9,5a18)') '#', '', 'c'
  do i=1,m
     write(1,'(i10,5e18.10)') i, y(i)
  enddo
  close(1)
end subroutine wts
