subroutine doe(idist,np,ndv,x,v)

implicit none

integer  i,j
integer  idist,np,ndv,iunit
real     x(np,ndv),v(np,ndv)
real     x_Ha(int(np/2),ndv),v_Ha(int(np/2),ndv)
real     x_FH(int(np/2),ndv),v_FH(int(np/2),ndv)
real     x_OI(int(np/4),ndv),v_OI(int(np/4),ndv)
real     vnorm

real, dimension(:), allocatable :: xmin,xmax

namelist/VARIABLES_RANGE/  xmin,xmax

!----------------------------------------------------------------------------

allocate(xmin(ndv),xmax(ndv))

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=VARIABLES_RANGE)
close(iunit)

vnorm = 1.0

if(idist.eq.1) then
  call D_HSS(x,v,vnorm,xmin,xmax,np,ndv)
endif

if(idist.eq.2) then
  call B_HSS(x,v,vnorm,xmin,xmax,np,ndv)
endif

if(idist.eq.3) then
  call D_HSS(x_Ha,v_Ha,vnorm,xmin,xmax,np/2,ndv)
  call B_HSS(x_FH,v_FH,vnorm,xmin,xmax,np/2,ndv)

  do i=1,int(np/2)
    do j=1,ndv
      x(i,j) = x_Ha(i,j)
      v(i,j) = v_Ha(i,j)
      x(int(i+np/2),j) = x_FH(i,j)
      v(int(i+np/2),j) = v_FH(i,j)
    end do
  end do
end if

if(idist.eq.4) then
  call D_HSS(x,v,vnorm,xmin,xmax,np,ndv)
  v = 0.
endif

if(idist.eq.5) then
  call B_HSS(x,v,vnorm,xmin,xmax,np,ndv)
  v = 0.
endif

if(idist.eq.6) then
  call D_HSS(x_Ha,v_Ha,vnorm,xmin,xmax,np/2,ndv)
  call B_HSS(x_FH,v_FH,vnorm,xmin,xmax,np/2,ndv)

  do i=1,int(np/2)
    do j=1,ndv
      x(i,j) = x_Ha(i,j)
      v(i,j) = v_Ha(i,j)
      x(int(i+np/2),j) = x_FH(i,j)
      v(int(i+np/2),j) = v_FH(i,j)
    end do
  end do
  v = 0.
end if

if(idist.eq.7) then
  call Ortho_Init(x_OI,v_OI,vnorm,xmin,xmax,np/4,ndv)

  do i=1,int(np/4)
    do j=1,ndv
      x(i,j)               = x_OI(i,j)*(-1.)
      v(int(i+np/4),j)     = v_OI(i,j)*(-1.)
      x(int(i+np/2),j)     = x_OI(i,j)
      v(int(i+np*3/4),j)   = v_OI(i,j)
    end do
  end do

  do i=1,np
    do j=1,ndv
      x(i,j)=((xmin(j)+xmax(j))+x(i,j)*(xmax(j)-xmin(j)))/2.
      v(i,j)=((xmin(j)+xmax(j))+v(i,j)*(xmax(j)-xmin(j)))/2.
    end do
  end do
end if

if(idist.eq.8) then
  open(unit=333,file='doe.inp',status='old',recl=10000)
   read(333,*) np
   do i=1,np
     read(333,*) x(i,:)
   end do
  close(333)
end if

deallocate(xmin,xmax)

!--------------------------------------------------------------------

return
end

