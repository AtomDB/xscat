      real function interp_mie_func(x, xderv)

      implicit none

      real x, xderv(5)
      integer DustType, MantleType, Drude, DustType_parm, Drude_parm
      real E_tmp, a_tmp, m_tmp, theta_tmp, x_tmp, rho
      common/fillers/Drude, DustType, MantleType,
     $     rho, theta_tmp, E_tmp, a_tmp, m_tmp, x_tmp

      real mie_func
      external mie_func
      integer i

      do i=1,5
         xderv(i) = 0.0e0
      end do
      DustType_parm = DustType
      Drude_parm = Drude
c      write(*,*) 'Call : ', rho,E_tmp,a_tmp,10**x
      interp_mie_func = mie_func(Drude_parm, DustType_parm, rho, E_tmp,
     $     a_tmp, 10.e0**x)

      return
      end
