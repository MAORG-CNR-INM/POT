! ----------------------------------------------------------------------
! Multi-objective synchronous deterministic particle swarm optimization
! ----------------------------------------------------------------------

subroutine MODPSO(ialgo,ndv,nfunc,nbox)

implicit none

include 'mpif.h'
include 'mpi.cmn'

integer    ndim,mdim,npdim
parameter (ndim=100,npdim=10*ndim,mdim=6000)

integer    i,j,m,t,r
integer    iter,ialgo,iflag,icoef,idist,iwall
integer    niter,nswarm,ndv,nfunc,nbox,nrow,nParetoGlobal,nParetoLimit
integer    npart,nfevalmax
integer    lmin(1)

integer    isolv,iobj,iunit

real       kpso,c1,c2
logical    bosol

integer    ktot,indice(npdim,mdim),weight(nfunc)
real       xtot(mdim,ndim),ftot(mdim,nfunc),fdum(ndim),fbox(ndim)
real       xmpi(mpidim,mpidim),fmpi(mpidim,mpidim)
real       fmpi2(mpidim,mpidim),box(mpidim,mpidim)

integer, dimension(:),     allocatable:: nParetoPersonal
real,    dimension(:),     allocatable:: xmin,xmax
real,    dimension(:,:),   allocatable:: x,v0,v,f,xglobal,xpersonal,ParetoGlobal,AggregFunct,fconstr
real,    dimension(:,:,:), allocatable:: ParetoPersonal
real,    dimension(:,:),   allocatable:: xaux,faux

character*80     filen

namelist/MODPSO_PARAMETERS/ npart,nfevalmax,icoef,idist,iwall,iflag
namelist/VARIABLES_RANGE/   xmin,xmax

!--------------------------------------------------------------------------------------------

allocate(xmin(ndv),xmax(ndv))

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=MODPSO_PARAMETERS)
close(iunit)

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=VARIABLES_RANGE)
close(iunit)

nswarm = npart
niter  = int(nfevalmax/nswarm)

allocate(x(nswarm,ndv),v(nswarm,ndv),v0(nswarm,ndv),f(nswarm,nfunc))
allocate(xglobal(nswarm,ndv),xpersonal(nswarm,ndv))
allocate(ParetoGlobal(nfevalmax,ndv+nfunc),AggregFunct(nswarm,niter))
allocate(ParetoPersonal(nswarm,niter,ndv+nfunc))
allocate(nParetoPersonal(nswarm))
allocate(fconstr(nswarm,nbox))

filen='sampl.XX.XX.XX.XX.plt'
write(filen(7:8),'(i2.2)') ialgo
write(filen(10:11),'(i2.2)') idist
write(filen(13:14),'(i2.2)') icoef
write(filen(16:17),'(i2.2)') iwall
open(500,file=filen,status='unknown',recl=10000)

filen='paret.XX.XX.XX.XX.plt'
write(filen(7:8),'(i2.2)') ialgo
write(filen(10:11),'(i2.2)') idist
write(filen(13:14),'(i2.2)') icoef
write(filen(16:17),'(i2.2)') iwall
open(501,file=filen,status='unknown',recl=10000)

! ---------------------------------------	SETTING NUMBER OF PARTICLES AND ITERATIONS

iter = 0

! ---------------------------------------	ALLOCATING ARRAYS	--------------------

do i=1,nfunc
  weight(i) = 0.50
end do

AggregFunct     = 0.
nParetoPersonal = 1

! ---------------------------------	INITIALIZATION PARTICLES POSITION AND VELOCITY	----

call doe(idist,nswarm,ndv,x,v0)

ParetoGlobal   = 0.
ParetoPersonal = 0.

nParetoGlobal  = 0
nParetoPersonal= 0

!---Set the coefficient 

if(icoef.eq.1) then     !Campana
  kpso = 1.000
  c1   = 0.400
  c2   = 1.300
end if
if(icoef.eq.2) then     !Diez
  kpso = 0.990
  c1   = 0.330
  c2   = 0.660
end if        
if(icoef.eq.3) then     !Eberhart-Clerc
  kpso = 0.729
  c1   = 2.050
  c2   = 2.050
end if        
if(icoef.eq.4) then     !Clerc
  kpso = 0.721
  c1   = 1.655
  c2   = 1.655
end if
if(icoef.eq.5) then     !Trelea
  kpso = 0.600
  c1   = 1.700
  c2   = 1.700
end if      
if(icoef.eq.6) then     !Carlisle-Dozier
  kpso = 0.729
  c1   = 2.300
  c2   = 1.800
end if              
if(icoef.eq.7) then     !Peri-Tinti
  kpso = 0.754
  c1   = 2.837
  c2   = 1.597        
end if

! ---------------------------------------	BEGIN ITERATION 	--------------------

!bosol=.true.

do while(iter*nswarm.lt.nfevalmax)
  iter = iter + 1
 
! --- Parallel

  do i=1,nswarm
    do j=1,ndv
      xmpi(i,j) = x(i,j)
    end do
  end do

  call parallel(nfunc,xmpi,nswarm,fmpi,fdum,ndim,iter,box)

  if(myrank.eq.0) then
    write(*,*)
    write(*,*) "----------------------------------------------------------------------------"
  end if

  do i=1,nswarm
    do j=1,nfunc
      fdum(j) = fmpi(i,j)
    end do
    do j=1,nbox
      fbox(j) = box(i,j)
    end do
    if(myrank.eq.0) then
      print *,"iter,part,f(N): ",iter,i,fdum(1:nfunc)
    end if
    do j=1,nfunc
      f(i,j) = fdum(j)
    end do
    do j=1,nbox
      fconstr(i,j) = fbox(j)
    end do
  end do

  if(myrank.eq.0) then
    write(*,*) "----------------------------------------------------------------------------"
    write(*,*)
  end if

! --- end parallel

  do i=1,nswarm
    write(500,*) iter,i,(x(i,j),j=1,ndv),(v0(i,j),j=1,ndv),(f(i,j),j=1,nfunc),fconstr(i,:)
  end do

! ---------------------------------------	COMPUTING FUNCTIONS	------------------

  do i=1,nswarm
    do j=1,ndv
      ParetoGlobal(nParetoGlobal+i,j) = x(i,j)
      ParetoPersonal(i,nParetoPersonal(i)+1,j) = x(i,j)
    end do

    do j=1,nfunc
      ParetoGlobal(nParetoGlobal+i,j+ndv) = f(i,j)
      ParetoPersonal(i,nParetoPersonal(i)+1,j+ndv) = f(i,j)
    end do
  end do

! ---------------------------------------	COMPUTING GLOBAL PARETO FRONT	--------

  nParetoGlobal = nParetoGlobal+nswarm
  nrow = nfevalmax !niter
  call Dominance(nParetoLimit,nParetoGlobal,nrow,ndv,nfunc,ParetoGlobal)

! ---------------------------------------	COMPUTING PERSONAL PARETO FRONT	------

  do i=1,nswarm
    nParetoPersonal(i)=nParetoPersonal(i)+1
    call Dominance(nParetoLimit,nParetoPersonal(i),niter,ndv,nfunc,ParetoPersonal(i,:,:))
  end do

! ---------------------------------------	COMPUTING ATTRACTORS	----------------

! ---------------------------------------	GLOBAL	------------------------------

  do i=1,nswarm
    call Closest(nParetoGlobal,ndv,ParetoGlobal(1:nParetoGlobal,1:ndv),x(i,1:ndv),xglobal(i,1:ndv))
  enddo

! ---------------------------------------	PERSONAL	----------------------------

  if (iflag==0) then
    do i=1,nswarm
      call Closest(nParetoPersonal(i),ndv,ParetoPersonal(i,1:nParetoPersonal(i),1:ndv),x(i,1:ndv),xpersonal(i,1:ndv))
    end do
  end if

  if (iflag==1) then
    AggregFunct(:,:) = 0.
      do i=1,nswarm
        do j=1,nParetoPersonal(i)
          do r=1,nfunc
            AggregFunct(i,j) = AggregFunct(i,j)+weight(r)*ParetoPersonal(i,j,ndv+r)
          end do
        end do
      end do

      do i=1,nswarm
        lmin = minloc(AggregFunct(i,1:nParetoPersonal(i)))
        xpersonal(i,:) = ParetoPersonal(i,lmin(1),1:ndv)
      end do
  end if

! ---------------------------------------	COMPUTING VELOCITY	------------------

  do i=1,nswarm
    do j=1,ndv
      v(i,j) = kpso*(v0(i,j)+c1*(xpersonal(i,j)-x(i,j))+c2*(xglobal(i,j)-x(i,j)))
      x(i,j) = x(i,j) + v(i,j)
    end do
  end do

! ---------------------------------------	PENALIZING	--------------------------			
  do i=1,nswarm
    do j=1,ndv
      if(x(i,j).lt.xmin(j)) then
        x(i,j) = xmin(j)
        if(iwall.eq.1) then
          v(i,j) = 0.
        else
          v(i,j) = -(v(i,j)/(kpso*(c1+c2)))
        end if
      end if
      if(x(i,j).gt.xmax(j)) then
        x(i,j) = xmax(j)
        if(iwall.eq.1) then
          v(i,j) = 0.
        else
          v(i,j) = -(v(i,j)/(kpso*(c1+c2)))
        end if
      end if
    end do
  end do    

  v0 = v

  write(501,*) nParetoGlobal,'nParetoGlobal'
  write(501,*) iter,'niter'
  write(501,*)
  do i=1,nParetoGlobal
    write(501,*) ParetoGlobal(i,:)
  end do
  write(501,*)
  write(501,*)
end do ! ---------------------------------------	THE END	----------------------

if(myrank.eq.0) then
  write(*,*) '---------------------------------------------------------------------------'
  write(*,*) nParetoGlobal,'nParetoGlobal'
  write(*,*) '---------------------------------------------------------------------------'
end if

open(555,file='pareto.out',status='unknown',recl=10000)
  do i=1,nParetoGlobal
    write(555,*) ParetoGlobal(i,:)
  end do
close(555)

allocate(xaux(nParetoGlobal,ndv))
allocate(faux(nParetoGlobal,nfunc))

open(20,file='pareto.out',status='old',recl=10000)
  do i=1,nParetoGlobal
    read(20,*) xaux(i,1:ndv),faux(i,1:nfunc)
  end do
close(20)

call sort(faux,xaux,nParetoGlobal,ndv,nfunc)

open(20,file='pareto.out',status='unknown',recl=10000)
  do i=1,nParetoGlobal
    write(20,*) xaux(i,1:ndv),faux(i,1:nfunc)
  end do
close(20)

deallocate(xaux,faux)

deallocate(xmin,xmax,x,v,v0,f)
deallocate(xglobal,xpersonal)
deallocate(ParetoGlobal,AggregFunct)
deallocate(ParetoPersonal,nParetoPersonal,fconstr)

! ------------------------------------------------------------------- 
return
end 




!----

subroutine sort(f,x,n,ndv,nfunc)

implicit none

integer i,j,n,ndv,nfunc
real    f(1:n,nfunc),x(1:n,ndv)
real    temp(nfunc),temp1(ndv)

!write(*,*) nfunc,f(1,:)

do i=1, n
  do j=n, i+1, -1
    if (f(j-1,1).lt.f(j,1)) then
      temp(:)=f(j-1,:)
      f(j-1,:)=f(j,:)
      f(j,:)=temp(:)

      temp1(:)=x(j-1,:)
      x(j-1,:)=x(j,:)
      x(j,:)=temp1(:)

    end if
  end do
end do


return
end
