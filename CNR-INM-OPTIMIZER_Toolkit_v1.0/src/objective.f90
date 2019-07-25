!--------------------------------------------------------------------
! Implement here your objective function
!
! A single objective and a multi-objective analytical problem 
! are provided as examples
!--------------------------------------------------------------------
subroutine objective_function(ndv,ndim,fdum)

implicit none

integer              iobj,iunit,ndv,ndim,j
real                 fdum(ndim),x(ndv)

namelist/OBJECTIVE/  iobj

!---------------------------------------------------------------------

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=OBJECTIVE)
close(iunit)

open(unit=20,file="variables.inp",status="unknown")
   do j=1,ndv
     read(20,*) x(j)
   end do
close(20)

if(iobj.eq.1) then
! --- Three-hump Camel (ndv=2); f_min(0,0) = 0
  fdum(1) = 2.*x(1)**2-1.05*x(1)**4+x(1)**6/6.+x(1)*x(2)+x(2)**2 
end if        

if(iobj.eq.2) then
! --- Three-hump Camel (ndv=2); f_min(0,0) = 0
  fdum(1) = 2.*x(1)**2-1.05*x(1)**4+x(1)**6/6.+x(1)*x(2)+x(2)**2
! --- Beale (ndv=2); min_f(3,0.5) = 0 
  fdum(2) = (1.5-x(1)+x(1)*x(2))**2 + (2.25-x(1)+x(1)*x(2)**2)**2 + (2.625-x(1)+x(1)*x(2)**3)**2
end if


return
end

