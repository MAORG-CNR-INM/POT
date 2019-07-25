subroutine linesearch(ndv,Gbst,xg,nfunc,step,fcount,iter,np,nfevalmax)

implicit none

include 'mpif.h'
include 'mpi.cmn'

integer          ndim,mdim,npdim
parameter       (ndim=100,npdim=10*ndim,mdim=6000)

integer          i,j,fcount,qcount,k,ndir
integer          iter,np,nfevalmax,iprob,ifun
integer          iunit,ialgo,iaux(1)
integer          ndv,nfunc,nbox
integer          nfr,nobj
integer          dir(2)
integer          jndv,LS_iter

real             xmpi(mpidim,mpidim),fmpi(mpidim,mpidim),box(mpidim,mpidim)
real             fdum(ndim),fbox(ndim)

real             Gbst,step!,fdum
real             xg(ndv),x_temp(ndv*2,ndv),xeval(ndv)
real             faux(2*ndv),fbst
real             gam,step0,eps,theta

logical          bofeasible,bosol

real, dimension(:), allocatable :: xmin,xmax
real, dimension(:,:),allocatable:: f,fconstr

namelist/VARIABLES_RANGE/         xmin,xmax
namelist/PROBLEM_PARAMETERS/      ndv,nfunc,nbox
namelist/LS_PARAMETERS/           gam,step0,qcount,eps,theta

! ---------------------------------------------------------------

allocate(xmin(ndv),xmax(ndv))

bofeasible = .true.
ndir   = 2
step   = 0.

fdum = 1.e+33

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=VARIABLES_RANGE)
close(iunit)

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=PROBLEM_PARAMETERS)
close(iunit)

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=LS_PARAMETERS)
close(iunit)

allocate(fconstr(2*ndv,nbox))
allocate(f(2*ndv,nfunc))

! --- The step is domain dipendent

step = step0*(maxval(xmax)-minval(xmin))

LS_iter = 0

do i=1,2*ndv
  x_temp(i,:) = xg(:)
end do

bosol=.true.

line: do while(step.gt.eps)

  LS_iter = LS_iter + 1

  do j=1,ndv
    x_temp(j,j)     = xg(j) + step
    x_temp(ndv+j,j) = xg(j) - step
  end do

! --- Wall type approach

  do i=1,2*ndv
    do j=1,ndv
      if(x_temp(i,j).gt.xmax(j)) then
        x_temp(i,j) = xmax(j)
      end if
      if(x_temp(i,j).lt.xmin(j)) then
        x_temp(i,j) = xmin(j)
      end if
    end do
  end do

!---Parallel

   do i=1,2*ndv
     do j=1,ndv
       xmpi(i,j) = x_temp(i,j)
     end do
   end do

   call parallel(nfunc,xmpi,2*ndv,fmpi,fdum,ndim,LS_iter,box,bosol)

   do i=1,2*ndv
     do j=1,nfunc
       fdum(j) = fmpi(i,j)
     end do
     do j=1,nbox
       fbox(j) = box(i,j)
     end do
     if(myrank.eq.0) then
       print *,"LS_iter,part,f: ",LS_iter,i,fdum(1)
     end if
     do j=1,nfunc
       f(i,j)   = fdum(j)
     end do
     do j=1,nbox
       fconstr(i,j) = fbox(j)
     end do
   end do

!---end parallel

   do i=1,2*ndv
     write(500,*) LS_iter,i,(x_temp(i,j),j=1,ndv),f(i,1),fconstr(i,:)
   end do

   fcount = fcount + 2*ndv

   faux(:) = f(:,1)

!---Update of global best

   if(minval(faux).lt.Gbst) then
     Gbst = minval(faux(:))
     iaux(:) = minloc(faux(:))
     xg(:) = x_temp(iaux(1),:)
     write(61,*) Gbst,'LINE'
     write(400,*) iter*np+fcount,xg(:),Gbst
     return
   else if((iter*np+fcount).ge.nfevalmax) then
     return 
   end if

  step = step*theta

end do line

deallocate(f,fconstr)
deallocate(xmin,xmax)

! ----
return
end
