! ------------------------------------------------------------------------
! Dolphin Pod Optimization
! ------------------------------------------------------------------------

subroutine DPO(ialgo,ndv,nfunc,nbox)

implicit none

include 'mpif.h'
include 'mpi.cmn'

integer          ndim,mdim,npdim
parameter       (ndim=100,npdim=10*ndim,mdim=6000)

integer          ndv,nfunc,nbox
integer          ndolp,nfevalmax,iline,icoef,idist,iwall
integer          iunit,nshoal
integer          i,j,k,iter,ialgo

real             xmpi(mpidim,mpidim),fmpi(mpidim,mpidim),box(mpidim,mpidim)
real             fdum(ndim),fbox(ndim),fbest(mdim)

real             d_tau,csi,k0,r0,a0,q,p,alp,b
real             fglob,ff_min,ff_max,eps,Rff
logical          bosol

real             xglob(ndv),Rx(ndv)

real,dimension(:),    allocatable:: aux,aux2,f,m
real,dimension(:,:),  allocatable:: delta,phi
real,dimension(:,:),  allocatable:: r,rbest
real,dimension(:,:,:),allocatable:: u,ubest
real,dimension(:),    allocatable:: xmin,xmax
real,dimension(:,:),  allocatable:: x,v,fconstr,xbest

real,dimension(:,:),  allocatable:: II,KK,GG,AA,VL,VR
integer                             INFO,LWORK,LWMAX
parameter                          (LWMAX = 100000)
real                                WORK(LWMAX)
real,dimension(:),    allocatable:: WR,WI,vect_tau

character*80     filen

namelist/DPO_PARAMETERS/  ndolp,nfevalmax,iline,icoef,idist,iwall
namelist/VARIABLES_RANGE/ xmin,xmax

!------------------------------------------------------------------------

allocate(xmin(ndv),xmax(ndv))

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=DPO_PARAMETERS)
close(iunit)

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=VARIABLES_RANGE)
close(iunit)

nshoal = ndolp

allocate(x(nshoal,ndv),v(nshoal,ndv),xbest(nshoal,ndv))
allocate(f(nshoal),fconstr(nshoal,nbox))
allocate(aux(nshoal),aux2(nshoal))
allocate(delta(nshoal,ndv),phi(nshoal,ndv))
allocate(r(nshoal,nshoal),rbest(nshoal,nshoal))
allocate(u(nshoal,nshoal,ndv),ubest(nshoal,nshoal,ndv))

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

! ----

iter    = 0

fglob   = -1.e+32
eps     = 1.e-12

xglob = 0.
phi   = 0.
delta = 0.
r     = 0.
rbest = 0.
u     = 0.
ubest = 0.
xbest = 0.
aux   = 0.
aux2  = 0.
fbest = -1.e+32

! -----------------------------------------------------

do k=1,ndv
  Rx(k) = abs(xmax(k)-xmin(k))
end do

ff_max = -1.e+32
ff_min =  1.e+32

!---Initialization first personal and global best

call doe(idist,nshoal,ndv,x,v)

do k=1,ndv
  x(:,k) = x(:,k)/Rx(k)
  v(:,k) = v(:,k)/Rx(k)
end do

do i=1,nshoal
  do j=1,ndv
    xbest(i,j) = x(i,j)
  end do
  fbest(i) = -1.e+32
end do

fglob = -1.e+32

!---Set the coefficient 

if (icoef.eq.1) then ! Serani and Diez (2015) Ndv<10
  csi = 1.0
  q   = 1.0
  p   = 8.0
  alp = 0.5
  b   = 0.0
end if
if (icoef.eq.2) then ! Serani and Diez (2015) Ndv>=10
  csi = 0.1 
  q   = 0.1 
  p   = 4.0
  alp = 0.5
  b   = 0.0
end if


k0 = 2*(1.-alp)*q/float(nshoal)
a0 = 2*alp*q/float(nshoal)
r0 = b*sqrt(float(nshoal))

allocate(II(nshoal,nshoal))
allocate(KK(nshoal,nshoal))
allocate(GG(nshoal,nshoal))
allocate(AA(2*nshoal,2*nshoal))
allocate(VR(2*nshoal,2*nshoal))
allocate(VL(2*nshoal,2*nshoal))
allocate(WR(2*nshoal))
allocate(WI(2*nshoal))
allocate(vect_tau(2*nshoal))

II = 0.
KK =-1.
GG = 0.
AA = 0.
VL = 0.
VR = 0.
WR = 0.
WI = 0.
WORK = 0.
vect_tau = 0.

do j=1,nshoal
  II(j,j) = 1.
  KK(j,j) = float(nshoal) -1.
end do

GG = csi*II
KK = -k0*KK

AA(       1:nshoal  ,nshoal+1:2*nshoal) = II(:,:)
AA(nshoal+1:2*nshoal,       1:nshoal  ) =-KK(:,:)
AA(nshoal+1:2*nshoal,nshoal+1:2*nshoal) =-GG(:,:)

LWORK = -1
call sgeev('N','N',2*nshoal,AA,2*nshoal,WR,WI,VL,2*nshoal,VR,2*nshoal,WORK,LWORK,INFO)
LWORK = MIN(LWMAX,INT(WORK(1)))
call sgeev('N','N',2*nshoal,AA,2*nshoal,WR,WI,VL,2*nshoal,VR,2*nshoal,WORK,LWORK,INFO)

do j=1,2*nshoal
  vect_tau(j) = 2.*abs(WR(j))/((WR(j))**2+(WI(j))**2)
end do

d_tau = minval(vect_tau(:))/p

!write(*,*) d_tau

deallocate(II,KK,GG,AA,VL,VR)
deallocate(WR,WI,vect_tau)

!---DPO iteration start

!bosol=.true.

do while(iter*nshoal.lt.nfevalmax)

  iter = iter + 1

!---Parallel

  do i=1,nshoal
    do j=1,ndv
      xmpi(i,j) = x(i,j)*Rx(j)
    end do
  end do

  call parallel(nfunc,xmpi,nshoal,fmpi,fdum,ndim,iter,box)

  do i=1,nshoal
    do j=1,nfunc
      fdum(j) = -fmpi(i,j)
    end do
    do j=1,nbox
      fbox(j) = -box(i,j)
    end do
    if(myrank.eq.0) then
      print *,"iter,part,f: ",iter,i,fdum(1)
    end if
    do j=1,nfunc
      f(i)   = fdum(j)
    end do
    do j=1,nbox
      fconstr(i,j) = fbox(j)
    end do
  end do

!---end parallel

!---Update of personal and global best

  do j=1,nshoal
    if(f(j).gt.fbest(j)) then
      fbest(j) = f(j)
      do k=1,ndv
        xbest(j,k) = x(j,k)
      end do
    end if

    if(f(j).gt.fglob) then  
      fglob = f(j)
      do k=1,ndv
        xglob(k) = x(j,k)
      end do
    end if
  end do

  do i=1,nshoal
    write(500,*) iter,i,x(i,:)*Rx(:),v(i,:)*Rx(:),-f(i),-fbest(i),-fconstr(i,:)
  end do

  write(501,*) iter,iter*nshoal,xglob(:)*Rx(:),-fglob

  ff_max = minval(f)
  ff_min = maxval(fbest)

  Rff = abs(ff_max-ff_min)
  if(Rff.eq.0.0) Rff = eps
  
! --- parameters

  do j=1,nshoal
    do i=1,nshoal
      r(i,j)  = 0.
      rbest(i,j) = 0.
      do k=1,ndv
        r(i,j)     = r(i,j)     + (x(i,k)    -x(j,k))**2
        rbest(i,j) = rbest(i,j) + (xbest(i,k)-x(j,k))**2
      end do
      r(i,j)     = sqrt(r(i,j))     + eps
      rbest(i,j) = sqrt(rbest(i,j)) + eps
      do k=1,ndv
        u(i,j,k)     = (x(i,k)     - x(j,k))/r(i,j)
        ubest(i,j,k) = (xbest(i,k) - x(j,k))/rbest(i,j)
      end do
    end do
  end do

! --- 

  do j=1,nshoal
    do k=1,ndv
      do i=1,nshoal
        aux2(i) = (r(i,j)-r0)*u(i,j,k)
        aux(i) = (2.*(fbest(i)-f(j))/Rff)*ubest(i,j,k)/(1.+sqrt(rbest(i,j)))
      end do
      delta(j,k) = k0*sum(aux2(:))
      phi(j,k)   = a0*sum(aux(:))
    end do

    do k=1,ndv

      v(j,k) = (1.-csi*d_tau)*v(j,k) + d_tau*(delta(j,k)+phi(j,k))
      x(j,k) = x(j,k) + v(j,k)*d_tau
 
! --- Wall-type approach not/elastic

      if(x(j,k).gt.xmax(k)/Rx(k)) then
        x(j,k) = xmax(k)/Rx(k)
        if(iwall.eq.1) then
          v(j,k) = 0.
        end if
        if(iwall.eq.2) then
          v(j,k) = -v(j,k)
        end if
      end if

      if(x(j,k).lt.xmin(k)/Rx(k)) then
        x(j,k) = xmin(k)/Rx(k)
        if(iwall.eq.1) then
          v(j,k) = 0.
        end if
        if(iwall.eq.2) then
          v(j,k) = -v(j,k)
        end if
      end if

    end do
  end do
!
end do

deallocate(x,v,xbest)
deallocate(f,fconstr)
deallocate(aux,aux2)
deallocate(delta,phi)
deallocate(r,rbest)
deallocate(u,ubest)

! ---------------------------------------------------------------
return
end
