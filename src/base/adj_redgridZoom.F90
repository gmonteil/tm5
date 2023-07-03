!### macro's ###################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "tm5.inc"
!
!###############################################################################

module adj_redgridZoom
  implicit none

  private
  public  :: adj_uni2red, adj_red2uni_em, adj_red2uni

  contains

  subroutine adj_uni2red(region)
   ! transforms data from uniform grid to reduced grid
   ! written by mike botchev, march-june 1999
   ! modified by Maarten Krol, dec 2002
   !
   !
   ! see comments red2uni
   ! 
   !
   !     Note that the weighting with m_uni moved to BEFORE the adjoint advection
   !

   use dims,        only: im,jm,lm
   use redgridZoom, only: nred, clustsize, jred, imredj
   use global_data, only: mass_dat
   use MeteoData  , only : m_dat
   use chem_param,  only: ntracet
   use ParTools  ,  only: ntracetloc
   implicit none
   !______________________________IO_______________________________________
   integer,intent(in)                      :: region
   !______________________________local________________________________________
   real,dimension(:,:,:,:),pointer         :: rm, rxm,rym,rzm
   real,dimension(:,:,:),  pointer         :: m

   integer i,ie,is,j,l,lrg,redfact,n
   real summ,sumrm, rmm
   !______________________________start________________________________________

   m => m_dat(region)%data
   rm => mass_dat(region)%rm_t
   rxm => mass_dat(region)%rxm_t
   rym => mass_dat(region)%rym_t
   rzm => mass_dat(region)%rzm_t

   do lrg=1,nred(region)
      redfact=clustsize(lrg,region)
      j = jred(lrg,region)
      do l=1,lm(region) 
         do i =  imredj(lrg,region),1,-1
            ! the i cell will be distributed within the is:ie array section
            is = (i-1)*redfact + 1  
            ie = i*redfact
            rmm = m(i,j,l)/(ie-is+1)
            m(is:ie,j,l) = rmm
            do n=1,ntracetloc
               rmm = rm(i,j,l,n) 
               rm(is:ie,j,l,n)= rmm
               rmm = rxm(i,j,l,n) 
               rxm(is:ie,j,l,n)= rmm
               ! rym and rzm are always distributed uniformly:
               rmm = rym(i,j,l,n) 
               rym(is:ie,j,l,n)= rmm
               rmm = rzm(i,j,l,n) 
               rzm(is:ie,j,l,n)= rmm
            enddo
         enddo
      enddo
   enddo

   nullify(m)
   nullify(rm)
   nullify(rxm)
   nullify(rym)
   nullify(rzm)

  end subroutine adj_uni2red

  subroutine adj_red2uni_em(region)
   ! transforms data from reduced grid back to uniform grid
   ! written by mike botchev, march-june 1999
   ! adjoint implementation: !ADJ MK: march 2003
   !     The reduced grid works as follows:
   !     o  m and am arrays are saved
   !     o  rm/rxm/rym/rzm arrays are compressed as a sum e.g.
   !         rm(2) = sum(rm(im/2+1:im))
   !         rm(1) = sum(rm(1:im)
   !     o  the array am is 'compressed', i.e. am(im/2+1) is moved to position 2
   !     o  advection is performed on the reduced im (2 in this example)
   !     o  rm/rxm/rym/rzm are reconstructed: 
   !         this used to be based on the advected m_uni, but because of the
   !         iteration in the advection, m is equally distributed (as rm, rxm,
   !         rym, rzm). Thus:
   !         mass = m(i,j,l) (with i the reduced grid counter)
   !         m(is:ie,j,l) = mass/(ie-is+1)  : put equal masses
   !         rm(im/2+1:im) = m(im/2+1:im)/mass   x   rm(2)
   !         rm(1:im/2+1)  = m(1:imi/2)/mass     x   rm(1)
   !         etc...
   !         Note that the forward run destroys the m information.
   !         m in the forward run remains on a 'x-mixed state' with equal masses
   !         over the reduced grid. This is important to know for the adjoint
   !         inplementation. After y- and z-advection we also mix the state of
   !         m, rm,rxm,rym,rzm acording to the reduced grid. 
   !         We have to treat this also in the adjoint code, in the advection
   !         routines itself. For now, we know that the boxes are ALWAYS mixed
   !         after an advection step in the forward run. Since other processes
   !         do not change m either, no problems with the adjoint are expected.
   !
   !   The adjoint implementation sould then be:
   !     o save am
   !     o adjoint red2uni_em:
   !         rm(1) = sum(rm(1:im/2))/(im/2)
   !         rm(2) = sum(rm(im/2+1:im)/(im/2)
   !     o reduce am
   !     o advection on reduced grid
   !     o adj_uni2red:
   !         rm(im/2+1:im) = rm(2)
   !         rm(1:im/2) = rm(1)
   !     o m should be restored by applying am and mixing in reduced grid. It
   !     should be checked if this is the case without doing something!
   !
   !
   use dims,        only: im,jm,lm
   use redgridZoom, only: nred, clustsize, jred, imredj
   use global_data, only: mass_dat
   use MeteoData  , only : m_dat
   use chem_param , only: ntracet
   use ParTools  ,  only: ntracetloc
   implicit none
   !______________________________IO_______________________________________
   integer,intent(in)                 :: region
   !______________________________local________________________________________
   real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
   real,dimension(:,:,:)  ,pointer    :: m
   integer i,ie,ii,is,j,l,lrg,n,redfact
   real summ,sumrm
   character*5 distr_mode
   !______________________________start________________________________________

   m => m_dat(region)%data
   rm => mass_dat(region)%rm_t
   rxm => mass_dat(region)%rxm_t
   rym => mass_dat(region)%rym_t
   rzm => mass_dat(region)%rzm_t

   do lrg=1,nred(region)
      redfact=clustsize(lrg,region)
      j = jred(lrg,region)

      do l=1,lm(region)
         do n = 1, ntracetloc
            ! JFM: add adjoint of periodic boundary condition 
            rm(im(region),j,l,n) = rm(im(region),j,l,n) + rm(0,j,l,n)
            rm(0,j,l,n) = 0.
            rxm(im(region),j,l,n) = rxm(im(region),j,l,n) + rxm(0,j,l,n)
            rxm(0,j,l,n) = 0.
            rym(im(region),j,l,n) = rym(im(region),j,l,n) + rym(0,j,l,n)
            rym(0,j,l,n) = 0.
            rzm(im(region),j,l,n) = rzm(im(region),j,l,n) + rzm(0,j,l,n)
            rzm(0,j,l,n) = 0.
         end do
         do i =  1,imredj(lrg,region)
            ! the is:ie  array section will be reduced to i
            is = (i-1)*redfact + 1  
            ie = i*redfact
            summ = sum(m(is:ie,j,l))
            m(i,j,l) = summ
            do n=1,ntracetloc
               sumrm = sum(rm(is:ie,j,l,n))/(ie-is+1)
               rm(i,j,l,n)  = sumrm
               sumrm = sum(rxm(is:ie,j,l,n))/(ie-is+1)
               rxm(i,j,l,n) = sumrm
               sumrm = sum(rym(is:ie,j,l,n))/(ie-is+1)
               rym(i,j,l,n) = sumrm
               sumrm = sum(rzm(is:ie,j,l,n))/(ie-is+1)
               rzm(i,j,l,n) = sumrm
            enddo   !n
         enddo  !i
         ! JFM: set remaining masses to zero
         do i = imredj(lrg,region)+1, im(region)
            do n = 1, ntracetloc
               rm(i,j,l,n) = 0.
               rxm(i,j,l,n) = 0.
               rym(i,j,l,n) = 0.
               rzm(i,j,l,n) = 0.
            end do  !n
         end do  !i
      enddo  !l
   enddo   !redgrid...

   nullify(m)
   nullify(rm)
   nullify(rxm)
   nullify(rym)
   nullify(rzm)

  end subroutine adj_red2uni_em

  subroutine adj_red2uni(region)
   ! transforms data from reduced grid back to uniform grid
   ! written by mike botchev, march-june 1999
   !
   ! JFM: follow old approach via m_uni
   !  o  need to advect air masses backward on uniform grid
   !     and forward on reduced grid to recover reduced-grid
   !     air masses as they were used in forward mode for
   !     distribution of tracer masses and slopes
   !
   use dims,        only: im,jm,lm
   use redgridZoom, only: nred, clustsize, jred, imredj, uni2red_mf, imred
   use global_data, only: mass_dat, wind_dat
   use MeteoData  , only : m_dat
   use chem_param , only: ntracet
   use ParTools  ,  only: ntracetloc
   implicit none
   !______________________________IO_______________________________________
   integer,intent(in)                 :: region
   !______________________________local________________________________________
   real,dimension(:,:,:,:),pointer    :: rm,rxm,rym,rzm
   real,dimension(:,:,:)  ,pointer    :: m,am
   real,dimension(-1:im(region)+2,-1:jm(region)+2,lm(region)) :: mb
   integer i,ie,ii,is,j,l,lrg,n,redfact,imr,jmr,lmr
   real summ,sumrm,sumrxm,sumrym,sumrzm
   character*5 distr_mode
   !______________________________start________________________________________

   m => m_dat(region)%data
   rm => mass_dat(region)%rm_t
   rxm => mass_dat(region)%rxm_t
   rym => mass_dat(region)%rym_t
   rzm => mass_dat(region)%rzm_t
   am => wind_dat(region)%am_t

   imr = im(region) ; jmr = jm(region) ; lmr = lm(region)

   ! First advect m backwards on uniform grid --> mb
   mb(1:imr,1:jmr,1:lmr)=m(1:imr,1:jmr,1:lmr) - am(0:imr-1,1:jmr,1:lmr) &
                                               +am(1:imr,  1:jmr,1:lmr)
   ! Reduce am
   call uni2red_mf(region)
   ! Coarsen mb to reduced grid
   do lrg=1,nred(region)
      redfact=clustsize(lrg,region)
      j = jred(lrg,region)
      do l=1,lmr
         do i =  1,imredj(lrg,region)
            ! the is:ie  array section will be reduced to i
            is = (i-1)*redfact + 1
            ie = i*redfact
            summ = sum(mb(is:ie,j,l))
            mb(i,j,l) = summ
         end do
      end do
   end do
   ! Advect forward on reduced grid
   do l = 1, lmr
      do j = 1, jmr
         do i = 1, imred(j,region)
            mb(i,j,l) = mb(i,j,l) + am(i-1,j,l) - am(i,j,l)
         end do
      end do
   end do

   ! Now do the actual adjoint of red2uni
   do lrg=1,nred(region)
      redfact=clustsize(lrg,region)
      j = jred(lrg,region)

      do l=1,lmr
         do n = 1, ntracetloc
            ! adjoint of periodic boundary condition 
            rm(imr,j,l,n) = rm(imr,j,l,n) + rm(0,j,l,n)
            rm(0,j,l,n) = 0.
            rxm(imr,j,l,n) = rxm(imr,j,l,n) + rxm(0,j,l,n)
            rxm(0,j,l,n) = 0.
            rym(imr,j,l,n) = rym(imr,j,l,n) + rym(0,j,l,n)
            rym(0,j,l,n) = 0.
            rzm(imr,j,l,n) = rzm(imr,j,l,n) + rzm(0,j,l,n)
            rzm(0,j,l,n) = 0.
         end do
         do i =  1,imredj(lrg,region)
            ! the is:ie  array section will be reduced to i
            is = (i-1)*redfact + 1  
            ie = i*redfact
            do n=1,ntracetloc
               sumrm = 0.
               sumrxm = 0.
               sumrym = 0.
               sumrzm = 0.
               do ii = is, ie
                  sumrm = sumrm + rm(ii,j,l,n)*m(ii,j,l)
                  sumrxm = sumrxm + rxm(ii,j,l,n)*m(ii,j,l)
                  sumrym = sumrym + rym(ii,j,l,n)*m(ii,j,l)
                  sumrzm = sumrzm + rzm(ii,j,l,n)*m(ii,j,l)
               end do
               rm(i,j,l,n)  = sumrm/mb(i,j,l)
               rxm(i,j,l,n)  = sumrxm/mb(i,j,l)
               rym(i,j,l,n)  = sumrym/mb(i,j,l)
               rzm(i,j,l,n)  = sumrzm/mb(i,j,l)
            enddo   !n
         enddo  !i
         ! JFM: set remaining masses to zero
         do i = imredj(lrg,region)+1, imr
            do n = 1, ntracetloc
               rm(i,j,l,n) = 0.
               rxm(i,j,l,n) = 0.
               rym(i,j,l,n) = 0.
               rzm(i,j,l,n) = 0.
            end do  !n
         end do  !i
      enddo  !l
   enddo   !redgrid...

   ! Finally copy mb to m
   m(1:imr,1:jmr,1:lmr) = mb(1:imr,1:jmr,1:lmr)

   nullify(m)
   nullify(rm)
   nullify(rxm)
   nullify(rym)
   nullify(rzm)
   nullify(am)

  end subroutine adj_red2uni

end module adj_redgridZoom
