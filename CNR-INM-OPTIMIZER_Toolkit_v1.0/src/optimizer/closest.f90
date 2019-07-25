Subroutine Closest(nPareto,ndv,Pareto,x,xClosest)

implicit none

integer:: ndv,nPareto,i,j,lmin(1)

real:: dmin

logical:: mask 

real,dimension(ndv):: x,xClosest
real,dimension(nPareto,ndv):: Pareto
real,dimension(nPareto,ndv+1):: Dist

Dist=01.e32

! --------------------------------------- COMPUTING DISTANCES	------------------

do i=1,nPareto
  Dist(i,ndv+1)=0.
  do j=1,ndv
    Dist(i,j)=Pareto(i,j)
    Dist(i,ndv+1)=Dist(i,ndv+1)+(Pareto(i,j)-x(j))**2
  end do
  Dist(i,ndv+1)=sqrt(Dist(i,ndv+1))
end do

! --------------------------------------- FINDING THE MINIMUM GLOBAL DISTANCE FOR EACH PARTICLE	

dmin = minval(Dist(:,ndv+1))
lmin = minloc(Dist(:,ndv+1))

xClosest=Dist(lmin(1),1:ndv)

!-----------------------------------------------------------------------

return
end
