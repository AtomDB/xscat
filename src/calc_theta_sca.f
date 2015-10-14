      real function calc_theta_sca(theta_obs, x)
c
c     Input is observed theta (in arcsec) and relative position x.
c     If small-angle formula is applicable, just use it.
c     Else, use exact value -- which is numerically unstable in large
c     parts of small-angle regions, as checked by a study with IDL.
c
c     Returns actual value in arcsec.
c
      implicit none

      include 'xscat.inc'

      real theta_obs, x

      real*8 std, ctd, ttd, xd, numer, denom

      calc_theta_sca = 0.0

      if ((theta_obs .lt. 300).or.
     $     (theta_obs .lt. 1000).and.(x.lt.0.99)) then ! Use small angle
         calc_theta_sca = theta_obs/(1.0-x)
      else
         std = sin(dble(theta_obs/radtoasec))
         ctd = cos(dble(theta_obs/radtoasec))
         ttd = tan(dble(theta_obs/radtoasec))
         xd  = dble(x)
c
c     RKS Solution
c         
c         numer = ctd*ctd - 2.0d0*xd + xd*xd/(ctd*ctd)
c         denom = xd*xd*ttd*ttd + (1.0d0-xd)*(1.0d0-xd)
c         calc_theta_sca = radtoasec*acos(sqrt(numer/denom))

c
c     Lynne Valencic Solution (more stable)
c         
         calc_theta_sca = radtoasec*atan(std/(ctd-xd))
c
c         write(*,*) 'Theta check: ',  theta_obs/(1.0-x),
c     $        radtoasec*acos(sqrt(numer/denom))
      end if
      
      return 
      end

      real function square(x, result, N)
      
      real x(*), result(*)
      integer i, N

      do i=1,N
         result(i) = x(i)*x(i)
      end do
      
      square = 0.0
      return
      end
