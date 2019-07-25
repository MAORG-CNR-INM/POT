! ----------------------------------------------------------------------
! Single-objective synchronous deterministic particle swarm optimization
! ----------------------------------------------------------------------
subroutine SODPSO(ialgo,ndv,nfunc,nbox)
! ----------------------------------------------------------------------

implicit none

include 'mpif.h'
include 'mpi.cmn'

integer          ndim,mdim,npdim
parameter       (ndim=100,npdim=10*ndim,mdim=6000)

integer          ndv,nfunc,nbox
integer          npart,nfevalmax,iline,icoef,idist,iwall
integer          iunit,nswarm
integer          i,j,k,iter,ialgo
integer          icount,qcount,fcount

real             xg(ndv)
real             xmpi(mpidim,mpidim),fmpi(mpidim,mpidim),box(mpidim,mpidim)
real             fdum(ndim),fbox(ndim),fbst(mdim,nfunc)
real             Gbst,chi,c1,c2
real             gam,step,step0,theta,eps
real             r_coef,lchar

real, dimension(:), allocatable :: xmin,xmax
real, dimension(:,:), allocatable :: x,v,f,fconstr,xbst

character*80     filen
character*4      success
character*1      morph

logical          update,update_LS,bosol

namelist/SODPSO_PARAMETERS/ npart,nfevalmax,iline,icoef,idist,iwall
namelist/VARIABLES_RANGE/   xmin,xmax
namelist/LS_PARAMETERS/     gam,step0,qcount,eps,theta

!----------------------------------------------------------------------------

allocate(xmin(ndv),xmax(ndv))

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=SODPSO_PARAMETERS)
close(iunit)

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=VARIABLES_RANGE)
close(iunit)

nswarm = npart

allocate(x(nswarm,ndv),v(nswarm,ndv),xbst(nswarm,ndv))
allocate(f(nswarm,nfunc))
allocate(fconstr(nswarm,nbox))


filen='sampl.XX.XX.XX.XX.XX.plt'
write(filen(7:8),'(i2.2)') ialgo
write(filen(10:11),'(i2.2)') idist
write(filen(13:14),'(i2.2)') icoef
write(filen(16:17),'(i2.2)') iwall
write(filen(19:20),'(i2.2)') iline
open(500,file=filen,status='unknown',recl=500)

filen='globa.XX.XX.XX.XX.XX.plt'
write(filen(7:8),'(i2.2)') ialgo
write(filen(10:11),'(i2.2)') idist
write(filen(13:14),'(i2.2)') icoef
write(filen(16:17),'(i2.2)') iwall
write(filen(19:20),'(i2.2)') iline
open(501,file=filen,status='unknown',recl=500)


!---LineSearch SET-UP

if(iline.eq.1) then
  open(iunit,file="POT.nml",action="read",recl=10000)
    read(iunit,NML=LS_PARAMETERS)
  close(iunit)
  open(unit=400,file="linegloba.plt",status="unknown",recl=10000)
end if

success = 'DPSO'
fcount  = 0
icount  = 0
update_LS = .true.

!---

iter = 0

!---Initialization first personal and global best

call doe(idist,nswarm,ndv,x,v)

do i=1,nswarm
  do j=1,ndv
    xbst(i,j) = x(i,j)
  end do
  fbst(i,1) = 1e+33
end do

Gbst = 1e+33

!---Set the coefficient 

if (icoef.eq.1) then ! Eberhart R, Shi Y, 2000; used in Diez M, Peri D, 2010 (single objective)
  chi = 0.729
  c1  = 2.050
  c2  = 2.050
endif
if (icoef.eq.2) then ! Carlisle-Dozier, 2001
  chi = 0.729
  c1  = 2.300
  c2  = 1.800
endif
if (icoef.eq.3) then ! Trelea, 2003
  chi = 0.600
  c1  = 1.700
  c2  = 1.700
endif
if (icoef.eq.4) then ! Clerc, 2007
  chi = 0.721
  c1  = 1.655
  c2  = 1.655
endif
if (icoef.eq.5) then ! Peri-Tinti, 2012
  chi = 0.754
  c1  = 2.837
  c2  = 1.597
endif

!---
!bosol=.true.

do while(iter*nswarm+fcount.lt.nfevalmax)
  iter = iter + 1
  update = .false.
  success= 'NONE'

!---Parallel

   do i=1,nswarm
     do j=1,ndv
       xmpi(i,j) = x(i,j)
     end do
   end do

   call parallel(nfunc,xmpi,nswarm,fmpi,fdum,ndim,iter,box)

   do i=1,nswarm
     do j=1,nfunc
       fdum(j) = fmpi(i,j)
     end do
     do j=1,nbox
       fbox(j) = box(i,j)
     end do
     if(myrank.eq.0) then
       print *,"iter,part,f: ",iter,i,fdum(1)
     end if
     do j=1,nfunc
       f(i,j)   = fdum(j)
     end do
     do j=1,nbox
       fconstr(i,j) = fbox(j)
     end do
   end do

!---end parallel

   do i=1,nswarm
     write(500,*) iter,i,(x(i,j),j=1,ndv),(v(i,j),j=1,ndv),f(i,1),fconstr(i,:)
   end do

!---Update of personal and global best

   do i=1,nswarm

!---Personal

     if(f(i,1).lt.fbst(i,1)) then
       fbst(i,1) = f(i,1)
       do j=1,ndv
         xbst(i,j) = x(i,j)
       end do
     end if
           
!---Global
     
     if(f(i,1).lt.Gbst) then  
       Gbst = f(i,1)
       update = .true.
       update_LS = .true.
       success= 'DPSO'
       do j=1,ndv
         xg(j) = x(i,j)
       end do
     end if
     write(61,*) Gbst,'DPSO'
   end do

! --- LineSearch after qcount without PSO improvement
 
   if(iline.eq.1) then
     if(.not.update) then
       if(update_LS) then
         icount = icount + 1
         if(icount.eq.qcount) then
           call linesearch(ndv,Gbst,xg,nfunc,step,fcount,iter,nswarm,nfevalmax)
           if(step.le.eps) then
             update_LS = .false.
             write(*,*) 'xg stationary point'
             success = 'NONE'
           else
             update_LS = .true.
             success = 'LINE'
           end if
           icount = 0
         end if
       end if
     end if
   end if

! ---

   write(501,*) iter,iter*nswarm+fcount,(xg(j),j=1,ndv),Gbst,success
!   success = 'DPSO'

!---DPSO

   do i=1,nswarm
     do j=1,ndv 
       v(i,j) = chi*(v(i,j)+c1*(xbst(i,j)-x(i,j))+c2*(xg(j)-x(i,j)))
       x(i,j) = x(i,j) + v(i,j)

! --- Wall type approach

       if(x(i,j).gt.xmax(j)) then
         x(i,j) = xmax(j)
         if(iwall.eq.1) then
           v(i,j) = 0.
         else
           v(i,j) = -v(i,j)/(chi*(c1+c2))
         end if
       end if
       if(x(i,j).lt.xmin(j)) then
         x(i,j) = xmin(j)
         if(iwall.eq.1) then
           v(i,j) = 0.
         else
           v(i,j) = -v(i,j)/(chi*(c1+c2))
         end if
       end if
        
     end do
   end do !end of DPSO
 end do  !end of while

write(502,*) iter,(xg(j),j=1,ndv),"fmin:",Gbst
close(61)
if(iline.eq.1) close(400)

deallocate(x,v,xbst,f,fconstr)
deallocate(xmax,xmin)

! ------------------------------------------------------------------- 
return
end
