      real function mie_interp(theta, result, N)
c
c     Uses the previously calculated Mie interpolation results
c     for our particular values of DustSize (a) and Energy (E)
c     to calculate the scattering cross section value for x,
c     assuming a fixed observing theta, which has an actual
c     scattering angle of theta/(1-x).
c
      implicit none

      include 'xscat.inc'

      integer N, i
      real theta,result, factor
      real interp_function, calc_theta_sca
      external interp_function, calc_theta_sca

      integer Drude, DustType, MantleType
      real theta0, distsq
      real rho, theta_tmp, E_tmp, a_tmp, m_tmp, x_tmp
      common/fillers/Drude, DustType, MantleType,
     $     rho, theta_tmp, E_tmp, a_tmp, m_tmp, x_tmp

      mie_interp = 0.0e0
      factor = 1.0
c
c If scattering angle < 0.001, just set to 0.001 to avoid numerical issues
c Scattering completely flat here anyway.
c
      if (theta .lt. minTheta) theta = minTheta
      result = interp_function(theta)

      if (result .ne. result) then ! check for NaNs
         write(*,*) 'Bad value at ', theta
         mie_interp = 1.0
      end if

c         write(*,*) 'func: ',theta0, x(i), factor, result(i)
      result = factor*result
c      write(*,*) theta, result

      return
      end
