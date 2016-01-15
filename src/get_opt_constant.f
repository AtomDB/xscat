      subroutine get_opt_constant(drude, dusttype, rho, energy, m)

      implicit none

      integer NMAXOPT,NDTYPES
      parameter(NMAXOPT=10000,NDTYPES=2)

      integer dusttype, drude, slot, ii
      real rho, energy, n_ana
      complex m

      complex interp_optconst
      external interp_optconst

      integer numopts(NDTYPES), dusttype_opt(NDTYPES)
      real e_opt(NMAXOPT,NDTYPES), n_opt(NMAXOPT,NDTYPES)
      real k_opt(NMAXOPT,NDTYPES)
      common/opt_parms/e_opt,n_opt,k_opt,numopts, dusttype_opt
      
      slot = -1
      m = cmplx(0.0,0.0)  ! default is vacuum

c     See Mauche & Gorenstein (1986, ApJ, 302, 372) for this definition.
      n_ana = 6.229e-4 * (rho/3.0) / energy**2

      if ((dusttype.eq.0).or.(drude.eq.1)) then  ! use Drude
         m = cmplx(n_ana,0.00)
      else
         do ii = 1, NDTYPES
            if (dusttype.eq.dusttype_opt(ii)) slot = ii
         end do
         if (slot .lt. 0) then 
            write(*,*) 'Dust type not loaded but should be.  Bad.'
            stop
         end if
         m = interp_optconst(slot, energy)
      end if

      return
      end

      subroutine init_opt_constant(slot, dusttype)

      implicit none

      integer NMAXOPT,NDTYPES
      parameter(NMAXOPT=10000,NDTYPES=2)
      real wksp(NMAXOPT)
      integer iwksp(NMAXOPT)
      integer dusttype, slot

      real e, n, k
      real e_in(NMAXOPT), n_in(NMAXOPT), k_in(NMAXOPT)

      integer numopts(NDTYPES), dusttype_opt(NDTYPES)
      real e_opt(NMAXOPT,NDTYPES), n_opt(NMAXOPT,NDTYPES)
      real k_opt(NMAXOPT,NDTYPES)
      common/opt_parms/e_opt,n_opt,k_opt,numopts, dusttype_opt
      
      character*160 filestem, filename, env
      character*80 dummy
      integer i, unit

c      if (numopts(slot) .gt. 0) return   ! Only load slot once

      call getenv('XSCAT', env)

      numopts(slot) = 0
      dusttype_opt(slot) = dusttype

      filestem='File not set'
      unit = 10
      if (dusttype .eq. 0)  filestem ='Drude Approximation'
      if (dusttype .eq. 1)  filestem ='/inputs/WD01_silicate_D03upd.dat'
      if (dusttype .eq. 2)  filestem ='/inputs/WD01_graphite_D03upd.dat'
      if (dusttype .eq. 3)  filestem ='/inputs/ZDA04_silicate.dat'
      if (dusttype .eq. 4)  filestem ='/inputs/ZDA04_amorphcar.dat'
      if (dusttype .eq. 5)  filestem ='/inputs/ZDA04_graphite.dat'
      if (dusttype .eq. 6)  filestem ='/inputs/ZDA04_refract.dat'
      if (dusttype .eq. 7)  filestem ='/inputs/ZDA04_CGS.dat' ! composite
      if (dusttype .eq. 8)  filestem ='/inputs/ZDA04_CGF.dat' ! composite
      if (dusttype .eq. 9)  filestem ='/inputs/ZDA04_CGB.dat' ! composite
      if (dusttype .eq. 10) filestem ='/inputs/ZDA04_CAS.dat' ! composite
      if (dusttype .eq. 11) filestem ='/inputs/ZDA04_CAF.dat' ! composite
      if (dusttype .eq. 12) filestem ='/inputs/ZDA04_CAB.dat' ! composite
      if (dusttype .eq. 13) filestem ='/inputs/ZDA04_CNS.dat' ! composite
      if (dusttype .eq. 14) filestem ='/inputs/ZDA04_CNF.dat' ! composite
      if (dusttype .eq. 15) filestem ='/inputs/ZDA04_CNB.dat' ! composite

      if (dusttype .eq. 100)
     $     filestem = '/inputs/vacuum.dat'
      if (dusttype .eq. 101)
     $     filestem = '/inputs/ice_h2o_optconst.dat'
      if (dusttype .eq. 102)
     $     filestem = '/inputs/ice_co2_optconst.dat'
      if (dusttype .eq. 103)
     $     filestem = '/inputs/ice_co_optconst.dat'

      if (dusttype .eq. 0) then 
c         write(*,FMT='("Using Drude Approximation (",I4,")")') dusttype
         numopts(slot) = 1
         return
      end if

      filename = trim(env)//filestem
c      write(*,*) ' '
c      write(*,FMT='("Using dust type",I4,":",A)') dusttype, filename

      open(unit=unit, err=100, FILE=filename, STATUS="OLD")
      if ((dusttype .ge. 3).and.(dusttype.le.6)) then ! ZDA orig
         do i=1,8 
            read(unit, *, end=1000, err=100) dummy
         end do 
      end if

      do i=1,NMAXOPT
         read(unit, *, end=1000, err=100) e, n, k
         numopts(slot) = numopts(slot) + 1

         if ((dusttype .ge. 3).and.(dusttype.le.6)) then  ! ZDA orig
            e_in(numopts(slot)) = 12.398e-4 / e  ! in keV
         else
            e_in(numopts(slot)) = e
         end if

         n_in(numopts(slot)) = n + 1.0

         if (k .gt. 0) k = -k   ! Use negative imaginary amplitudes
         k_in(numopts(slot)) = k
      end do
      
 100  continue
      write(*,*) 'Error reading optical constant file:',filename,
     $     'for dusttype =', dusttype
      stop

 1000 continue
      close(unit)

      call sort3(numopts, e_in, n_in, k_in, wksp, iwksp)
         
      do i=1,numopts(slot)
         e_opt(i, slot) = e_in(i)
         n_opt(i, slot) = n_in(i)
         k_opt(i, slot) = k_in(i)
      end do

      return
      end


      complex function interp_optconst(slot, energy)
C
      implicit none

      integer NMAXOPT,NDTYPES
      parameter(NMAXOPT=10000,NDTYPES=2)

      integer slot
      real energy

      integer jl, ju, jm
      real GRAD,DF

      integer numopts(NDTYPES), dusttype_opt(NDTYPES)
      real e_opt(NMAXOPT,NDTYPES), n_opt(NMAXOPT,NDTYPES)
      real k_opt(NMAXOPT,NDTYPES)
      common/opt_parms/e_opt,n_opt,k_opt,numopts, dusttype_opt

      real n_res, k_res
      complex result

      result = cmplx(1.0,0.0) ! default is vacuum

      if ((energy.lt.e_opt(1,slot)).or.
     $     (energy.gt.e_opt(numopts(slot),slot))) then
         write(*,*) 'Energy = ',energy,
     $        ' beyond range of optical constants for dustype ',
     $        dusttype_opt(slot)
         interp_optconst = result
         return 
      end if
C     
      JL=0
      JU=numopts(slot)+1

 10   if (JU-JL.gt.1) then
         JM=(JU+JL)/2
         if (energy .gt. e_opt(JM,slot)) then 
            JL=JM
         else 
            JU=JM
         end if
         go to 10 
      end if 
C     
C------energy is sandwiched between JL and JU ------
C
      if(n_opt(JL,slot).gt.0.0.and.n_opt(JU,slot).gt.0.) then
         GRAD=(LOG10(n_opt(JU,slot))-LOG10(n_opt(JL,slot)))/
     $        (LOG10(e_opt(JU,slot))-LOG10(e_opt(JL,slot)))
         DF=GRAD*(LOG10(energy)-LOG10(e_opt(JL,slot)))
         n_res = 10**(LOG10(n_opt(JL,slot))+DF)
      else
         n_res = n_opt(JL,slot)+
     $        (n_opt(JU,slot)-n_opt(JL,slot))*
     $        ((energy-e_opt(JL,slot))/(e_opt(JU,slot)-e_opt(JL,slot)))
      end if

      if(k_opt(JL,slot).gt.0.0.and.k_opt(JU,slot).gt.0.) then
         GRAD=(LOG10(k_opt(JU,slot))-LOG10(k_opt(JL,slot)))/
     $        (LOG10(e_opt(JU,slot))-LOG10(e_opt(JL,slot)))
         DF=GRAD*(LOG10(energy)-LOG10(e_opt(JL,slot)))
         k_res = 10**(LOG10(k_opt(JL,slot))+DF)
      else
         k_res = k_opt(JL,slot)+
     $        (k_opt(JU,slot)-k_opt(JL,slot))*
     $        ((energy-e_opt(JL,slot))/(e_opt(JU,slot)-e_opt(JL,slot)))
      end if

      result = cmplx(n_res, k_res)
      interp_optconst = result

      return
      end
      
