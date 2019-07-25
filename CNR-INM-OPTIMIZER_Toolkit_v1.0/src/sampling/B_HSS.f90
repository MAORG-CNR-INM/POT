      subroutine B_HSS(xp,v,vnorm,xdvmin,xdvmax,n,ndv)

! -----------------------------------------------------------------------

      parameter(nprimi=2000)
      real p(nprimi)
      real a(0:nprimi)
      real x(1000)
      real xdvmin(ndv),xdvmax(ndv)
      real xp(n,ndv)
      real v(n,ndv)
      real ddum,xcorr,vnorm

      external generaprimi

! ----------------------------------------------------------------------

      call prime(p,17390)

      ndvp1   = ndv + 1
      if(ndvp1.gt.2000) write(*,*) 'Low number of prime numbers'

! -- nprimiU va tarato!
      nprimiU = 20
        
      do iface=1,ndv*2

        do i=1,ceiling(1.*n/ndv/2.0) !indice di punti per faccia

          x(1) = float(i-1)/float(n-1)*(2*ndv)
          jaux = 1

          do j=2,ndvp1

            jbase = j

            xcorr = float(i)

            do k=nprimiU,1,-1

              ddum = p(jbase)**k
              if(ddum.le.xcorr) then
                idum = int(xcorr/ddum)
                a(k) = float(idum)
                xcorr = xcorr - a(k)*ddum
              else
                a(k) = 0.d0
              end if

            end do

            a(0) = xcorr

            x(j) = 0.
            do k=1,nprimiU
              ddum = p(jbase)**k
              x(j) = x(j) + (a(k)/ddum)
            end do

            ipoint=((iface-1)*n/ndv/2)+i
            jpoint=j-1
            if (iface.eq.(2*jpoint-1)) then
              xp(ipoint,jpoint)= 0.
              v(ipoint,jpoint)=0.
            else if (iface.eq.(2*jpoint)) then
              xp(ipoint,jpoint)= 1.
              v(ipoint,jpoint) = 1.     
            else
              xp(ipoint,jpoint)=x(jaux)
              v(ipoint,jpoint) = -(x(jaux)-0.5)*2.0/sqrt(float(ndv))
              jaux=jaux+1
            endif
          end do
        end do
      end do

      do i=1,n
        do j=1,ndv
           xp(i,j) = xdvmin(j) + xp(i,j)*(xdvmax(j)-xdvmin(j))
           v(i,j) =  xdvmin(j) + v(i,j)*(xdvmax(j)-xdvmin(j))
        end do
      end do

! ----------------------------------------------------------------------
      return
      end
! ----------------------------------------------------------------------
