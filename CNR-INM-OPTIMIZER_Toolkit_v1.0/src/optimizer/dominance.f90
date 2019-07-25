subroutine Dominance(nParetoLimit,nPareto,nrow,ndv,nfunct,Pareto)

implicit none

integer:: nrow,i,j,nPareto,niter,ncount,ifunct,nParetoTemp,ndv,nfunct,nParetoLimit
logical:: bopareto

real,dimension(nrow,ndv+nfunct):: ParetoTemp,Pareto
real,dimension(ndv):: x,xClosest
real,dimension(nPareto,nPareto):: Dist

!----------------------------------------------------------------------------------

ParetoTemp = 0.
ParetoTemp = Pareto
Pareto     = 0.

nParetoTemp=nPareto

nPareto = 0

Dist=0.

do i = 1,nParetoTemp
  bopareto = .true.
  j = 1
  do while((bopareto).and.(j.le.nParetoTemp))
    ncount = 0
    do ifunct = 1,nfunct 
      if (ParetoTemp(j,ndv+ifunct).le.ParetoTemp(i,ndv+ifunct)) then
         ncount = ncount+1
      endif
    enddo
    if(ncount==nfunct) then
      do ifunct = 1,nfunct 
        if (ParetoTemp(j,ndv+ifunct).lt.ParetoTemp(i,ndv+ifunct)) then
          ncount = ncount+1
        endif
      enddo
    endif
    if(ncount.gt.nfunct) then
      bopareto = .false.
    end if
    j = j+1
  enddo
  if(bopareto) then
    nPareto = nPareto+1
    Pareto(nPareto,:) = ParetoTemp(i,:)
  endif
enddo

ParetoTemp = 0.
ParetoTemp = Pareto
Pareto     = 0.

nParetoTemp= nPareto
nPareto    = 1
Pareto(1,:)= ParetoTemp(1,:)

outer: do i = 2,nParetoTemp
         do j = 1,nPareto
           if(Pareto(j,ndv+1) == ParetoTemp(i,ndv+1)) then
! -- Found a match so start looking again
             cycle outer
           end if
         end do
! -- No match found so add it to the output
         nPareto = nPareto + 1
         Pareto(nPareto,:) = ParetoTemp(i,:) 
       end do outer

!------------------------------------------------------------------------------

return
end

