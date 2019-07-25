Program main

implicit none

include 'mpif.h'
include 'mpi.cmn'

real t_start,t_end,t_exec

character*6 dire
character*12 command1
character*18 command2

integer iunit,ierrmpi,clock
integer ndv,nfunc,nbox
integer ialgo

namelist/PROBLEM_PARAMETERS/      ndv,nfunc,nbox
namelist/OPTIMIZATION_PARAMETERS/ ialgo

! **********************************************************************

call system_clock(count=clock)
t_start = clock

! ----------------------------------------------------------------------
! --  MPI startup
! ----------------------------------------------------------------------

call mpi_init(ierrmpi)

! -- Set nprocs if serial run

call mpi_comm_size(mpi_comm_world,nprocs,ierrmpi)
call mpi_comm_rank(mpi_comm_world,myrank,ierrmpi)

if(nprocs.eq.0) nprocs = 1

if(myrank.eq.0) then
  write(*,*)
  write(*,*) '------------------------------------------------------------------------- '
  write(*,*) '|                                                                        |'
  write(*,*) '|                        Parallel Optimizer Toolkit                      |'
  write(*,*) '|                                   v1.0                                 |'
  write(*,*) '|                       24 Jul. 2019 ... release 1.0                     |'
  write(*,*) '|   Multidisciplinary Analysis and Optimization Research Group (MAORG)   |'
  write(*,*) '|                           @CNR-INM, Rome, Italy                        |'
  write(*,*) '|             Serani A., Pellegrini R., Leotardi C., Diez M.             |'
  write(*,*) '|                                                                        |'
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)
end if

! ----------------------------------------------------------------------
! -- Initialize Dirs
! ----------------------------------------------------------------------

write(*,*) 'Number of processors: ',nprocs,myrank
command1 = 'mkdir CPU000'
write(command1(10:12),'(i3.3)') myrank
call system(command1)

command2 = 'cp -fr *.* CPU000/'
write(command2(15:17),'(i3.3)') myrank
call system(command2)

dire = 'CPU000'
write(dire(4:6),'(i3.3)') myrank
call chdir(dire)

! ----------------------------------------------------------------------

iunit = 20

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=PROBLEM_PARAMETERS)
close(iunit)

open(iunit,file="POT.nml",action="read",recl=10000)
  read(iunit,NML=OPTIMIZATION_PARAMETERS)
close(iunit)
call optimizer(ialgo,ndv,nfunc,nbox)

! ----------------------------------------------------------------------
! -- MPI closure
! ----------------------------------------------------------------------

call mpi_finalize(ierrmpi)

call system_clock(count=clock)
t_end = clock
t_exec = (t_end-t_start)/10000

write(33,*) "Elapsed CPU time = ",t_exec,"seconds"

end program


