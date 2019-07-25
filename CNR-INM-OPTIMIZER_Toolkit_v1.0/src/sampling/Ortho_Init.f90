subroutine Ortho_Init(x,v,vnorm,xdvmin,xdvmax,nswarm,ndv)

implicit none

integer nswarm,ndv
integer i,j

real xdvmin(ndv),xdvmax(ndv)
real x(nswarm,ndv)
real v(nswarm,ndv)
real vnorm

do i=1,nswarm
  x(i,i) = 1.
  v(i,i) = 1.
end do

return
end 
