!------------------------------------------------------------------------------
subroutine optimizer(ialgo,ndv,nfunc,nbox)
!------------------------------------------------------------------------------

implicit none

include 'mpif.h'
include 'mpi.cmn'

integer ialgo,ndv,nfunc,nbox

!-----------------------------------------------------------------------


if(ialgo.eq.1) then
  call SODPSO(ialgo,ndv,nfunc,nbox)
end if

if(ialgo.eq.2) then
  call DPO(ialgo,ndv,nfunc,nbox)
end if

if(ialgo.eq.3) then
  call MODPSO(ialgo,ndv,nfunc,nbox)
end if

!-----------------------------------------------------------------------
return
end

