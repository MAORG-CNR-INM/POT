c ----------------------------------------------------------------------
      subroutine parallel(nfunc,xmpi,npts,fmpi,fdum,ndim,iter,box)
c ----------------------------------------------------------------------

      include 'mpif.h'
      include 'mpi.cmn'

      integer           istart(mpidim),iend(mpidim)
      integer           master,npts,ierrmpi
      integer           iunit
      integer           ndv,nfunc,nbox

      real              xmpi(mpidim,mpidim),fmpi(mpidim,mpidim)
      real              box(mpidim,mpidim)
      real              xdum(mpidim,mpidim),fdum(ndim)

      logical           bofeasible,bosol

      character*6       dire

      namelist/PROBLEM_PARAMETERS/      ndv,nfunc,nbox

c-----------------------------------------------------------------------

      bofeasible = .true.
      iunit = 20

      open(iunit,file="POT.nml",action="read",recl=10000)
        read(iunit,NML=PROBLEM_PARAMETERS)
      close(iunit)

c ----------------------------------------------------------------------

      if(myrank.eq.0) then
        open(unit=22,file='computandi.sunt',
     &       status='unknown',form='formatted',recl=500)
        do i=1,npts
          write(22,200) (xmpi(i,j),j=1,ndv)
        end do
        close(22) 
      end if

      call mpi_barrier(mpi_comm_world,ierrmpi)

c-- set i-start and i-end for parallel

      call set_isviev(istart,iend,npts)
      call mpi_barrier(mpi_comm_world,ierrmpi)

c -- objective functions evaluation

      write(*,*)
      write(*,*) '-----------------------------------------------------'

      do 100 i=istart(myrank+1),iend(myrank+1)

c -- myrank is the processor number (0:n-1)

        write(*,*) 'Solution n. ',i

        open(unit=20,file="variables.inp",status="unknown")
        do j=1,ndv
          write(20,*) xmpi(i,j)
        end do
        close(20)

        do j=1,mpidim
          fmpi(i,j)=0.
        end do

        call objective_function(ndv,ndim,fdum)        

        do j=1,ndim
          fmpi(i,j) = fdum(j)
        end do

 100  continue

      call mpi_barrier(mpi_comm_world,ierrmpi)

c -- updating rank = 0 (from rank > 0)

      call ofun_par(xdum,istart,iend,xmpi,mpidim,mpidim)
      call mpi_barrier(mpi_comm_world,ierrmpi)

      call ofun_par(xdum,istart,iend,fmpi,mpidim,mpidim)
      call mpi_barrier(mpi_comm_world,ierrmpi)

      call ofun_par(xdum,istart,iend,box,mpidim,mpidim)
      call mpi_barrier(mpi_comm_world,ierrmpi)


c-- NEW BROADCAST

      master = 0

      if(nprocs.gt.1) then
        do j=1,mpidim
          call mpi_bcast(fmpi(1,j),mpidim,MPI_REAL,master,
     &                   MPI_COMM_WORLD,ierrmpi)
          call mpi_barrier(mpi_comm_world,ierrmpi)
        end do
        do j=1,mpidim
          call mpi_bcast(box(1,j),mpidim,MPI_REAL,master,
     &                   MPI_COMM_WORLD,ierrmpi)
          call mpi_barrier(mpi_comm_world,ierrmpi)
        end do
        do j=1,ndv
          call mpi_bcast(xmpi(1,j),mpidim,MPI_REAL,0,
     &                   MPI_COMM_WORLD,ierrmpi)
          call mpi_barrier(mpi_comm_world,ierrmpi)
        end do
      end if

c --

  200 format(1000e15.7)

c-----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine ofun_par(xdum,isv,iev,x,iidim,jjdim)

c-----------------------------------------------------------------------

      include 'mpif.h'
      include 'mpi.cmn'

      integer isv(mpidim),iev(mpidim)
      real x(iidim,jjdim)
      real xdum(iidim,jjdim)

c-----------------------------------------------------------------------

      if(nprocs.gt.1) then
        if(myrank.ge.1) call ofun_sen(xdum,isv,iev,x,iidim,jjdim)
        if(myrank.eq.0) call ofun_rec(xdum,isv,iev,x,iidim,jjdim)
      end if

c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------

      subroutine ofun_sen(sen_pack,isv,iev,x,iidim,jjdim)

c-----------------------------------------------------------------------

      include 'mpif.h'
      include 'mpi.cmn'

      integer isp,icount,idest
      integer isv(mpidim),iev(mpidim)
      real sen_pack(iidim,jjdim),x(iidim,jjdim)

c-----------------------------------------------------------------------

      isp=1
      do i=isv(myrank+1),iev(myrank+1)
        do j=1,jjdim
          sen_pack(isp,j)=x(i,j)
        end do
        isp=isp+1
      end do
      icount=iidim*jjdim
      idest=0
      call mpi_send(sen_pack,icount,mpi_real,idest,nprocs*myrank, 
     &              mpi_comm_world,ierrmpi)

      call mpi_barrier(mpi_comm_world,ierrmpi)

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------

      subroutine ofun_rec(rec_pack,isv,iev,x,iidim,jjdim)

c-----------------------------------------------------------------------

      include 'mpif.h'
      include 'mpi.cmn'

      integer icount,source,isp,idest
      integer istat(mpi_status_size),isv(mpidim),iev(mpidim)
      real rec_pack(iidim,jjdim),x(iidim,jjdim)

c-----------------------------------------------------------------------

      icount=iidim*jjdim
      do source=1,nprocs-1
        call mpi_recv(rec_pack,icount,mpi_real,source,nprocs*source,
     &                mpi_comm_world,istat,ierrmpi)
        isp=1
        do i=isv(source+1),iev(source+1)
          do j=1,jjdim
            x(i,j)=rec_pack(isp,j)
          end do
          isp=isp+1
        end do
      end do

      call mpi_barrier(mpi_comm_world,ierrmpi)

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------

      subroutine set_isviev(isv,iev,lens)

c-----------------------------------------------------------------------

      include 'mpif.h'
      include 'mpi.cmn'

      integer lens,disie
      integer ske 
      integer isv(mpidim),iev(mpidim)

c-----------------------------------------------------------------------

      disie=lens/nprocs
      ske=lens-disie*nprocs

c-- for balanced loads

      do i=1,nprocs
        isv(i)=disie*(i-1)+1
        iev(i)=isv(i)+disie-1
      end do

c-- for unbalanced loads

      if(ske.gt.0) then
        do i=1,ske
          isv(i)=isv(i)+i-1
          iev(i)=isv(i)+disie
        end do
        do i=ske+1,nprocs
          isv(i)=isv(i)+ske
          iev(i)=iev(i)+ske
        end do
      end if

      if(myrank.eq.0) then
        open(99,file='smista.out',status='unknown',recl=500)
        write(99,1001)
        do i=1,nprocs
          write(99,1002)i-1,isv(i),iev(i),iev(i)-isv(i)+1
        end do
        close(99)
      end if

1001  format(t2,'rank',t11,'i-start',t22,'i-end',
     &       t32,'delta is-ie') 
1002  format(t2,i2,t11,i5,t22,i5,
     &       t32,i5)

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
