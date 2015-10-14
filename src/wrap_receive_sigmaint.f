      subroutine wrap_receive_sigmaint(Energy, size_a, size_m, xpos,
     $     Interp, Norm, Drude_in, DustType_in, MantleType_in,
     $     rho_in, Tmax, epsilon, result)
      
c
c     Goal is to calculate 
c     \sigma_{LOS}(E, \phi)  = \int_{{\phi}\over{1-x}}^{\pi} 
c              2 \pi \theta {{d\sigma}\over{d\Omega}}(E, \theta, a, x) d\theta
c     
c
      implicit none

      include 'xscat.inc'

      integer ScatType, Interp, Drude_in
      integer DustType_in, MantleType_in
      real Energy, size_a, size_m, rho_in, Norm, Xmin, Xmax
      real xpos
      real Tmin, Tmax
      real epsilon
      real ActualTsca, FullTsca
      integer iAng, numT

      real Intensity, result
 
      real mie_func
      external mie_func

c      real coatmie_func
c      external coatmie_func

      real interp_mie_func, mie_interp
      external interp_mie_func, mie_interp

c      real interp_coatmie_func
c      external interp_coatmie_func

c      real mie_x, mie_mx
c      external mie_x, mie_mx

      real calc_theta_sca
      external calc_theta_sca

      integer Drude, DustType, MantleType
      real rho, theta_tmp, E_tmp, a_tmp, m_tmp, x_tmp
      common/fillers/Drude, DustType, MantleType,
     $     rho, theta_tmp, E_tmp, a_tmp, m_tmp, x_tmp

      Drude = Drude_in
      DustType = DustType_in 
      MantleType = MantleType_in
      rho = rho_in
      E_tmp = Energy
      a_tmp = size_a
      m_tmp = size_m

c
c     Initialize the optical constants
c
      call init_opt_constant(1, DustType)
      if (size_m .gt. 0.0) call init_opt_constant(2, MantleType)
c
c     Now we need to figure out the type of run.  Is there 
c     a mantle on these grains, or are they homogeneous?
c     If size_m > 0.0, then there is a mantle and we need to
c     use the coated mie scheme, rather than MIEV0.
c
      
      if (Interp.eq.1) then ! do interpolation
c         write(*,'(A)',advance='no') 'Started building approximation...'
         if (size_m .eq. 0.0) then ! no mantle
            call build_interp_function(interp_mie_func,
     $           minTheta, radtoasec, epsilon)
         else                      ! with mantle
            write(*,*) 'Do not use mantles yet.'
c            call build_interp_function(interp_coatmie_func,
c     $           minTheta, radtoasec, epsilon)
         end if
c         write(*,*) 'done.'
      else
         write(*,*) 'Use the interpolation.'
      end if

c
c
      ActualTsca = calc_theta_sca(Tmax, xpos)

      call qtrap(mie_interp, ActualTsca, radtoasec, Intensity, epsilon)

c      call qtrap(mie_interp, 0, ActualTsca, Intensity, epsilon)
c      call qtrap(mie_interp, 0.0, radtoasec, Intensity, epsilon)
c      write(*,*) 'Intensity = ', Intensity

      result = Intensity
      
      return
      end
