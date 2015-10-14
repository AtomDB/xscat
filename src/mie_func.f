      real function mie_func(drude, dusttype, rho, energy,
     $     dustradius, theta)

      implicit none

      integer maxang, momdim
      parameter(maxang=40, momdim=10)

      real XX, xmu(maxang)
      complex crefin
      integer drude, dusttype
      real energy     ! in keV
      real dustradius ! in um
      real theta      ! in arcseconds
      real intensity  ! in cm^2 sr^-1

      logical anyang, perfct, prnt(2)
      integer ipolzn, numang, nmom
      real gqsc, mimcut, pmom(0:momdim, 1000), qext, qsca, spike
      complex  sforw, sback, s1(maxang),s2(maxang), tforw(2), tback(2)
      
      real rho
      real radtoasec,pi, hc
c      external miev0

      perfct = .false.
      mimcut = 0.0e0
      anyang = .true.
      nmom = 0
      ipolzn = 0
      prnt(1) = .false.
      prnt(2) = .false.
      pi = 3.1415926535
      hc = 12.39843e-4 ! h*c, in units of keV um
      radtoasec = (180./pi)*60.*60.   ! 1 radian = 206265 arcsec

      call get_opt_constant(drude, dusttype, rho, energy, crefin)

      XX = 2*pi*dustradius/(hc/energy)

      numang = 1
      xmu(1) = cos(theta/radtoasec)

      call miev0( XX, crefin, perfct, mimcut, anyang,
     $     numang, xmu, nmom, ipolzn, momdim, prnt, qext, qsca,
     $     gqsc, pmom, sforw, sback, s1, s2, tforw, tback, spike)
c
c Intensity has units of cm^2 sr^-1
c

      intensity=0.5e0*(hc*1.e-4/(2.0*pi*Energy))**2
     $        *(abs(s1(1))**2 + abs(s2(1))**2)

c      write(*,*) theta, XX, Energy, dustradius, intensity
c      write(*,*) theta, intensity
      mie_func = intensity/radtoasec ! integration done in arcsec, not radians
      mie_func = 2*pi*(theta/radtoasec)*mie_func ! multiply by 2pi radius

      return
      end

