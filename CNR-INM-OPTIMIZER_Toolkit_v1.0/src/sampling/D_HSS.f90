      subroutine D_HSS(xp,v,vnorm,xdvmin,xdvmax,n,ndv)

!-----------------------------------------------------------------------

      parameter(nprimi=2000)
      real p(nprimi)
      real a(0:nprimi)
      real xcorr
      real x(7000)
      real xdvmin(ndv),xdvmax(ndv)
      real xp(n,ndv)
      real v(n,ndv),vnorm

      external generaprimi

! ----------------------------------------------------------------------

      call prime(p,17390)

      ndvp1   = ndv + 1
      if(ndvp1.gt.2000) write(*,*) 'Low number of prime numbers'

! -- nprimiU va tarato!

      nprimiU = 10

      do i=1,n

        x(1) = float(i-1)/float(n-1)

        do j=2,ndvp1

          jbase = j + 0

          xcorr = float(i-1)

          do k=nprimiU,1,-1

            if(p(jbase).le.xcorr) then
              idum = int(xcorr/(p(jbase)**k))
              a(k) = float(idum)
              xcorr = xcorr - a(k)*(p(jbase)**k)
            else
              a(k) = 0.
            end if

          end do

          a(0) = xcorr

          x(j) = 0.
          do k=0,nprimiU
            x(j) = x(j) + a(k)/(p(jbase)**(k+1))
          end do

        end do

        do j=1,ndv
!          xp(i,j) = x(j)
          xp(i,j) = xdvmin(j) + x(j)*(xdvmax(j)-xdvmin(j))
!          v(i,j)  = vnorm*(xp(i,j)+0.5)*2.0/sqrt(2.*ndv)
          v(i,j) = (xp(i,j)-(xdvmin(j)+xdvmax(j))/2.)*(2./sqrt(float(ndv)))
       end do

      end do

! ----------------------------------------------------------------------
      return
      end
! ----------------------------------------------------------------------
