
      SUBROUTINE MieLay( RCORE, RSHELL, WVNO, RINDSH, RINDCO, MU,
     &                   NUMANG, QEXT, Qsca, QBS, GQSC,
     &                   M1, M2, S21, D21, MaxAng )

c *******************************************************************
c    This subroutine computes electromagnetic scattering by a
c    stratified sphere (a particle with a spherical core surrounded
c    by a spherical shell).  The surrounding medium is assumed to
c    have refractive index unity.  The formulas, manipulated to avoid
c    the ill-conditioning that plagued earlier formulations, were 
c    published in:

c        Toon, O. and T. Ackerman, Applied Optics 20, 3657 (1981)

c    The latest version of this program is available by anonymous ftp 
c    from climate.gsfc.nasa.gov in directory pub/wiscombe.

c    The program was based on the famous homogeneous sphere
c    program of Dave, published in:

c       Dave, J.V., "Subroutines for Computing the Parameters of the
c          Electromagnetic Radiation Scattered by a Sphere",
c          IBM Scientific Center, Palo Alto, California,
c          Report No. 320 - 3236, May 1968

c       Dave, J.V., Applied Optics 8, 155 (1969)

c    This was done because the formulas are identical in structure to
c    those for the homogeneous sphere, except that the coefficients of
c    the Mie series (commonly denoted by little-a-sub-n and 
c    little-b-sub-n) are much more complicated.

c    The address of the first author is:

c         Dr. O. B. Toon (toon@sky.arc.nasa.gov)
c         NASA Ames Research Center
c         M.S. 245-3
c         Moffett Field, CA (USA)

c    The program was explicitly designed to avoid the ill-conditioning
c    which was latent in the standard analytic formulation of the
c    Mie core-shell problem.  Oddly, this ill-conditioning had been
c    exposed and eliminated in the late 1960s for the homogeneous sphere
c    problem, but went unrecognized for the core-shell problem, leading
c    to many incorrect results being published prior to 1981.  In
c    particular, previous calculations were generally wrong for the
c    case of a thin carbon shell.

c    After a number of years of experience, this program seems to have
c    only two limitations:  a slow degradation as size increases,
c    and a rapid degradation as size decreases due to lack of explicit
c    handling of the ill-conditioning in that limit. It has been used 
c    successfully for cases with large imaginary refractive index (both
c    in core and shell) and even with real refractive index less than 
c    unity.

c    For too-large particles, internal array sizes will be inadequate,
c    but this generates an error message and an abort.  

c    It is highly recommended to use the DOUBLE PRECISION version of 
c    this program, called 'DMiLay', on machines with 32-bit floating-
c    point arithmetic (e.g. all IEEE arithmetic machines).  'DMiLay' 
c    is also available on the network.  'MieLay' may be adequate but
c    this should be tested by running it back to back with 'DMiLay'.

c        Note that in 32-bit arithmetic it was impossible to run 
c        'MieLay' above shell radius = 3 (with WVNO=1, so shell size
c        parameter = 3 also) due to overflow, whereas 'DMiLay' can be
c        run well above this limit due to a larger range of exponents.

c    The original version of this program defaulted to a homogeneous
c    sphere case when the core was less than 10**(-6) the size of the
c    shell.  It could also have done so when core and shell radii were
c    equal although it did not.  But this option was dropped since 
c    it could better be tested for in the calling program;  if a
c    homogeneous sphere case is detected, it is far better to call one 
c    of the programs designed explicitly for that case.

c    NOTE:  This program requires input scattering angles between
c           zero and 90 degrees.  Then it does calculations for
c           those angles, plus all their supplements.  Thus, to get,
c           e.g., 170 degrees, you must use 10 degrees.

c    The program was modified and further documented by W. Wiscombe 
c    (wiscombe@climate.gsfc.nasa.gov; NASA Goddard, Code 913,
c    Greenbelt, MD 20771), including:
c    ** complex refractive indices of shell and core submitted as 2
c          complex arguments rather 4 real arguments
c    ** scattering angles submitted as cosines rather than angles
c          themselves (anticipating important usage where cosines remain
c          constant while program is called repeatedly for different
c          shell and/or core sizes in order to integrate over size)
c    ** returning scattering matrix elements M1, M2, D21, S21 separately
c          rather than bundled into awkward data structure ELTRMX
c    ** defining new arrays S1 and S2, the complex scattering 
c          amplitudes, for which ELTRMX had formerly been used as
c          temporary storage (allowing easy modification to return S1
c          and S2 through argument list rather than M1, M2, D21, S21)
c    ** use of internal work arrays, with guards against blowing their
c          dimensions, rather than submitting these arrays as arguments;
c          allows easy conversion to ALLOCATABLE arrays in Fortran-90
c    ** elimination of all GO TOs but one
C    ** elimination of dangerous EQUIVALENCE statement (it was used to
c          equate variables, not to share storage)
c    ** more mnemonic variable names
c    ** error-checking input arguments
C    ** making data structures more optimal for vectorization
c          (particularly PI, TAU, and the replacements for ELTRMX);
c          mainly so that innermost loops are over the first
c          dimension of all arrays
c    ** polishing and declaration standardization using NAG Fortran
c          Tool nag_decs
c    ** creation of a DOUBLE PRECISION version using NAG Fortran
c          Tool nag_apt
c    ** certification by 'flint' (Fortran 'lint', for C folks)

c    Suggestions for the future:  
c    ** much of this program reflects the belief, true in the late
c         1960s but not any longer, that complex arithmetic should
c         be broken into real and imaginary parts for efficiency,
c         in spite of massive loss of clarity;  the program would
c         simplify considerably if done in complex variables
c         (now only partially true)
c    ** improve treatment of angles greater than 90 degrees
c         (an awkward and confusing hangover from the old Dave code)
c    ** phrase program in terms of size parameters and eliminate
c         input argument WVNO
c    ** develop special-case formulas for Rcore-->0 for any Rshell,
c         and for Rshell-->0 (comparing single and double precision 
c         versions showed larger and larger differences as size
c         parameter fell below about 1.E-2, esp. for gQsc); the
c         layered sphere formulae are ill-conditioned in this limit 
c         just as the homogeneous sphere formulae are


c    I N P U T   A R G U M E N T S

c    (Definition:  size parameter = sphere circumference / wavelength )

c      Rshell      radius of shell

c      Rcore       radius of core

c      WVNO        2*pi / wavelength

c      RindSh      COMPLEX refractive index of shell (negative
c                     imaginary part)

c      RindCo      COMPLEX refractive index of core (negative
c                     imaginary part)

c      MU          array of cosines of scattering angles (angles between
c                     directions of incident and scattered radiation).
c                     For angles between 90 and 180 degrees, use the
c                     supplement (180-angle) of the angle instead, so
c                     0.le.MU.le.1 (see comments below on M1,M2,21,D21)

c      NumAng      Number of scattering angles for which computations
c                     are required; should not exceed MaxAng
C                     (NOTE:  NumAng=0 will suppress the calculation
c                      of the scattering matrix quantitities  M1, M2,
c                      S21, D21 and save a lot of computer time)

c      MaxAng      First dimension of M1,M2,21,D21 in calling program



c    O U T P U T   A R G U M E N T S

c      (Definitions for these arguments can be found in H.C. van de
c       Hulst, Light Scattering By Small Particles, Dover Press, New
c       York, 1981 (reprint of 1957 edition); abbreviated VDH below)

c      QEXT     Efficiency factor for extinction (VDH Sec 9.32)
c               (same as corresponding quantity in MIEV)

c      Qsca     Efficiency factor for scattering (VDH Sec 9.32)
c               (same as corresponding quantity in MIEV)

c      GQSC     average(cosine theta) * Qsca (VDH Sec 9.32)
c                  (<cos theta> is usually denoted by g, hence
c                   the name of the variable)
c               (same as corresponding quantity in MIEV)

c      QBS      Backscatter cross section.
c               ( Re(SBACK)**2 + Im(SBACK)**2 ) / (Pi*XSHELL**2)
c               where the corresponding quantity from MIEV is

c               SBACK = 0.5*sum(n=1 to inf)((-1)**(n+1)(2n+1)(an-bn))

c               and an,bn are ACOE,BCOE below.

c      M1(j,k)  Element M1 of scattering matrix F' (VDH Sec 5.14);
c                  M1(j,1) refers to angle with cosine MU(j); 
c                  M1(j,2) refers to supplement of that angle.
C               (Be sure to type REAL in calling program.)

c      M2(j,k)  Element M2 of scattering matrix F' (VDH Sec 5.14);
c                  M2(j,1) refers to angle with cosine MU(j); 
c                  M2(j,2) refers to supplement of that angle.
C               (Be sure to type REAL in calling program.)

c     S21(j,k)  Element S21 of scattering matrix F' (VDH Sec 5.14);
c                  S21(j,1) refers to angle with cosine MU(j); 
c                  S21(j,2) refers to supplement of that angle.

c     D21(j,k)  Element D21 of scattering matrix F' (VDH Sec 5.14);
c                  D21(j,1) refers to angle with cosine MU(j); 
c                  D21(j,2) refers to supplement of that angle.


c    L O C A L   V A R I A B L E S

c      ACOE     (COMPLEX) Mie coeff. little-a-sub-n
c      BCOE     (COMPLEX) Mie coeff. little-b-sub-n
c      ACOEM1   (COMPLEX) Mie coeff. little-a-sub-(n-1)
c      BCOEM1   (COMPLEX) Mie coeff. little-b-sub-(n-1)

c      LL       dimension of ACAP; in original program was 7000; for
c               conserving memory this should not be much bigger than
c                  1.1*Abs(RindSh) * x + 1
c               BUT it must ALWAYS exceed 150

c      TA(1,2)  real, imaginary parts of WFN(1)
c      TA(3,4)  real, imaginary parts of WFN(2)

c      S1,S2   complex scattering amplitudes; these are processed to
c               produce elements of the real scattering matrix
c               (they would need to be multiplied by a factor to
c               give the true scattering amplitudes, VDH Sec 9.31)

C      TOLER   tolerance for cosines of angles when slightly below 0
c               or slightly above 1 (a common occurrence)

c    NOTE:  the definitions of U(i) in this program are not the same as
c           the u-sub-i defined in the Toon/Ackerman paper.  The
c           correspondence is:
c             usub1 = u(1)    usub2 = u(5)
c             usub3 = u(7)    usub4 = dumsq
c             usub5 = u(2)    usub6 = u(3)
c             usub7 = u(6)    usub8 = u(4)
c             ratio of spherical Bessel to spherical Hankel func = u(8)

c    The Bessel function ratio A is always computed by downward 
c    recurrence.

c **********************************************************************

c     .. Parameters ..

      INTEGER   MxAng, LL
      REAL      ZERO, ONE, TWO
      PARAMETER ( MxAng = 100, LL = 1000, ZERO = 0.0, ONE = 1.0,
     &            TWO = 2.0 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   MaxAng, NUMANG
      REAL      GQSC, QBS, QEXT, Qsca, RCORE, RSHELL, WVNO
      COMPLEX   RINDCO, RINDSH
c     ..
c     .. Array Arguments ..

      REAL      MU( NUMANG ), D21( MAXANG, 2 ), M1( MAXANG, 2 ),
     &          M2( MAXANG, 2 ), S21( MAXANG, 2 )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR, PASS1
      INTEGER   I, J, K, M, N, NMX1, NMX2, NN

      REAL      AA, ARE, AIM, AM1RE, AM1IM, BB, BRE, BIM, BM1RE,
     &          BM1IM,CC, COSX1, COSX4, DD, DENOM, E2Y1, EY1, 
     &          EY1MY4, EY1PY4, EY4, PINUM, RMM, RX, SINX1, SINX4,
     &          TOLER, X1, X4, XCORE, XSHELL, Y1, Y4

      COMPLEX   AC, ACOE, ACOEM1, BC, BCOE, BCOEM1, CI, CZERO, 
     &          DH1, DH2, DH4, DUMMY, DUMSQ, K1, K2, K3,
     &          P24H21, P24H24, RRFX, SBACK, WM1
c     ..
c     .. Local Arrays ..

      REAL      PI( MxAng, 3 ), SI2THT( MxAng ), TAU( MxAng, 3 ),
     &          T( 4 ), TA( 4 )

      COMPLEX   ACAP( LL ), S1( MxAng, 2 ), S2( MxAng, 2 ), 
     &          U( 8 ), W( 3, LL ), WFN( 2 ), Z( 4 )
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, ASIN, CMPLX, COS, EXP, MOD, REAL, SIN
c     ..
c     .. Save statement ..

      SAVE  PINUM, PASS1

c     .. Data statements ..

      DATA      PASS1 / .True. / ,  TOLER / 1.E-6 /,
     &          CZERO / ( 0., 0. ) / ,  CI / ( 0., 1. ) /
c     ..


      IF( PASS1 ) THEN

         PINUM  = TWO*ASIN( ONE )
         PASS1  = .False.

      END IF


      XSHELL = RSHELL*WVNO
      XCORE  = RCORE*WVNO
      T( 1 ) = XSHELL*ABS( RINDSH )
      NMX1   = 1.10*T( 1 )
      NMX2   = T( 1 )
      IF( NMX1.LE.150 ) THEN
         NMX1  = 150
         NMX2  = 135
      END IF

c                        ** Check input arguments for gross errors

      INPERR = .False.

      IF( WVNO.LE.0.0 ) INPERR = WRTBAD( 'WVNO' )

      IF( RSHELL.LE.0.0 ) INPERR = WRTBAD( 'Rshell' )

      IF( RCORE.LE.0.0 .OR. RCORE.GT.RSHELL ) 
     &    INPERR = WRTBAD( 'Rcore' )

      IF( REAL(RINDSH).LE.0.0 .OR. AIMAG(RINDSH).GT.0.0 )
     &    INPERR = WRTBAD( 'RindSh' )

      IF( REAL(RINDCO).LE.0.0 .OR. AIMAG(RINDCO).GT.0.0 ) 
     &    INPERR = WRTBAD( 'RindCo' )

      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NumAng' )

      IF( NUMANG.GT.MxAng ) INPERR = WRTDIM( 'MxAng', NUMANG )

      IF( NUMANG.GT.MaxAng ) INPERR = WRTDIM( 'MaxAng', NUMANG )

      IF( NMX1+1 .GT. LL ) INPERR = WRTDIM( 'LL', NMX1 + 1 )

      DO 10 J = 1, NUMANG
         IF( MU(J).LT. -TOLER .OR. MU(J).GT.1.0+TOLER ) 
     &       INPERR = WRTBAD( 'MU' )
   10 CONTINUE

      IF( INPERR ) CALL ERRMSG(
     &             'MIELAY--Input argument errors.  Aborting...',.True.)


      K1     = RINDCO*WVNO
      K2     = RINDSH*WVNO
      K3     = CMPLX( WVNO )
      Z( 1 ) = RINDSH*XSHELL
      Z( 2 ) = XSHELL
      Z( 3 ) = RINDCO*XCORE
      Z( 4 ) = RINDSH*XCORE
      X1     =  REAL( Z(1) )
      Y1     = AIMAG( Z(1) )
      X4     =  REAL( Z(4) )
      Y4     = AIMAG( Z(4) )
      RX     = ONE / XSHELL

c                                ** Down-recurrence for A function
      ACAP( NMX1 + 1 ) = CZERO
      DO 20 M = 1, 3
         W( M, NMX1 + 1 ) = CZERO
   20 CONTINUE

      RRFX  = ONE / ( RINDSH*XSHELL)
      DO 40  NN = NMX1, 1, -1

         ACAP( NN ) = (( NN + 1 )*RRFX) -
     &                ONE / ( ((NN + 1)*RRFX) + ACAP(NN + 1) )

         DO 30 M = 1, 3

            W( M, NN ) = (( NN+1 ) / Z( M+1 )) -
     &                   ONE / ( (( NN+1 ) / Z(M+1)) + W( M, NN + 1) )
   30    CONTINUE

   40 CONTINUE


      DO 50 J = 1, NUMANG

         SI2THT( J ) = ONE - MU( J )**2
         PI( J, 1 )  = ZERO
         PI( J, 2 )  = ONE
         TAU( J, 1 ) = ZERO
         TAU( J, 2 ) = MU( J )

   50 CONTINUE

c                          ** Initialization of homogeneous sphere

      T( 1 ) = COS( XSHELL )
      T( 2 ) = SIN( XSHELL )
      WM1      = CMPLX( T(1), - T(2) )
      WFN( 1 ) = CMPLX( T(2), T(1) )
      TA( 1 ) = T( 2 )
      TA( 2 ) = T( 1 )
      WFN( 2 ) = RX*WFN( 1 ) - WM1
      TA( 3 ) =  REAL( WFN(2) )
      TA( 4 ) = AIMAG( WFN(2) )

c                      ** Initialization procedure for stratified sphere
      N      = 1
      SINX1  = SIN( X1 )
      SINX4  = SIN( X4 )
      COSX1  = COS( X1 )
      COSX4  = COS( X4 )
      EY1    = EXP( Y1 )
      E2Y1   = EY1**2
      EY4    = EXP( Y4 )
      EY1MY4 = EXP( Y1 - Y4 )
      EY1PY4 = EY1*EY4
      AA     = SINX4*( EY1PY4 + EY1MY4 )
      BB     = COSX4*( EY1PY4 - EY1MY4 )
      CC     = SINX1*( E2Y1 + ONE )
      DD     = COSX1*( E2Y1 - ONE )
      DENOM  = ONE + E2Y1*( 4.0*SINX1**2 - TWO + E2Y1 )
      DUMMY  = CMPLX( ( AA*CC + BB*DD ) / DENOM, 
     &                ( BB*CC - AA*DD ) / DENOM )
      DUMMY  = DUMMY *( ACAP(N) + N / Z(1) ) / ( W(3, N) + N / Z(4) )
      DUMSQ  = DUMMY**2

      P24H24 = 0.5 + CMPLX( SINX4**2 - 0.5, COSX4*SINX4 )*EY4**2
      P24H21 = 0.5*CMPLX( SINX1*SINX4 - COSX1*COSX4, 
     &                    SINX1*COSX4 + COSX1*SINX4 ) *EY1PY4
     &       + 0.5*CMPLX( SINX1*SINX4 + COSX1*COSX4, 
     &                  - SINX1*COSX4 + COSX1*SINX4 ) *EY1MY4
      DH1    = Z( 1 ) / ( ONE + CI*Z(1) ) - ONE / Z( 1 )
      DH2    = Z( 2 ) / ( ONE + CI*Z(2) ) - ONE / Z( 2 )
      DH4    = Z( 4 ) / ( ONE + CI*Z(4) ) - ONE / Z( 4 )
      P24H24 = P24H24 / (( DH4 + N / Z(4) )*( W(3, N) + N / Z(4) ))
      P24H21 = P24H21 / (( DH1 + N / Z(1) )*( W(3, N) + N / Z(4) ))

      U( 1 ) = K3*ACAP( N ) - K2*W( 1, N )
      U( 2 ) = K3*ACAP( N ) - K2*DH2
      U( 3 ) = K2*ACAP( N ) - K3*W( 1, N )
      U( 4 ) = K2*ACAP( N ) - K3*DH2
      U( 5 ) = K1*W( 3, N ) - K2*W( 2, N )
      U( 6 ) = K2*W( 3, N ) - K1*W( 2, N )
      U( 7 ) = - CI*( DUMMY*P24H21 - P24H24 )
      U( 8 ) = TA( 3 ) / WFN( 2 )

      ACOE  = U( 8 )*( U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5) ) /
     &               ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )

      BCOE  = U( 8 )*( U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6) ) /
     &               ( U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6) )

      ACOEM1  = ACOE
      BCOEM1  = BCOE
      ARE =  REAL( ACOE )
      AIM = AIMAG( ACOE )
      BRE =  REAL( BCOE )
      BIM = AIMAG( BCOE )

      QEXT  = 3.*( ARE + BRE )
      Qsca  = 3.*( ARE**2 + AIM**2 + BRE**2 + BIM**2 )
      GQSC  = ZERO
      SBACK = 3.*( ACOE - BCOE )
      RMM   = ONE

      AC  = 1.5*ACOE
      BC  = 1.5*BCOE
      DO 60 J = 1, NUMANG

         S1( J, 1 ) = AC*PI( J, 2 ) + BC*TAU( J, 2 )
         S1( J, 2 ) = AC*PI( J, 2 ) - BC*TAU( J, 2 )
         S2( J, 1 ) = BC*PI( J, 2 ) + AC*TAU( J, 2 )
         S2( J, 2 ) = BC*PI( J, 2 ) - AC*TAU( J, 2 )

   60 CONTINUE


c ***************** Start of Mie summing loop ******************

      N  = 2
   70 CONTINUE
c                           ** Recurrences for functions little-pi,
c                              little-tau of Mie theory
      T( 1 ) = 2*N - 1
      T( 2 ) = N - 1
      DO 80 J = 1, NUMANG

         PI( J, 3 ) = ( T(1)*PI(J, 2)*MU( J ) - N*PI(J, 1) ) / T(2)

         TAU( J, 3 ) = MU( J )*( PI(J, 3) - PI(J, 1) ) -
     &                 T( 1 )*SI2THT( J )*PI( J, 2 ) + TAU( J, 1 )

   80 CONTINUE

c                                 ** Here set up homogeneous sphere
      WM1      = WFN( 1 )
      WFN( 1 ) = WFN( 2 )
      WFN( 2 ) = T( 1 )*RX*WFN( 1 ) - WM1
      TA( 1 ) =  REAL( WFN(1) )
      TA( 2 ) = AIMAG( WFN(1) )
      TA( 3 ) =  REAL( WFN(2) )
      TA( 4 ) = AIMAG( WFN(2) )

c                                 ** Here set up stratified sphere

      DH1    = - N / Z( 1 ) + ONE / ( N / Z(1) - DH1 )
      DH2    = - N / Z( 2 ) + ONE / ( N / Z(2) - DH2 )
      DH4    = - N / Z( 4 ) + ONE / ( N / Z(4) - DH4 )
      P24H24 = P24H24 / (( DH4 + N / Z(4) )*( W(3, N) + N / Z(4) ))
      P24H21 = P24H21 / (( DH1 + N / Z(1) )*( W(3, N) + N / Z(4) ))
      DUMMY  = DUMMY *( ACAP(N) + N / Z(1) ) / ( W(3, N) + N / Z(4) )
      DUMSQ  = DUMMY**2

      U( 1 ) = K3*ACAP( N ) - K2*W( 1, N )
      U( 2 ) = K3*ACAP( N ) - K2*DH2
      U( 3 ) = K2*ACAP( N ) - K3*W( 1, N )
      U( 4 ) = K2*ACAP( N ) - K3*DH2
      U( 5 ) = K1*W( 3, N ) - K2*W( 2, N )
      U( 6 ) = K2*W( 3, N ) - K1*W( 2, N )
      U( 7 ) = - CI*( DUMMY*P24H21 - P24H24 )
      U( 8 ) = TA( 3 ) / WFN( 2 )

      ACOE  = U( 8 )*( U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5) ) /
     &               ( U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5) )

      BCOE  = U( 8 )*( U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6) ) /
     &               ( U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6) )
      ARE =  REAL( ACOE )
      AIM = AIMAG( ACOE )
      BRE =  REAL( BCOE )
      BIM = AIMAG( BCOE )

c                           ** Increment sums for efficiency factors

      AM1RE =  REAL( ACOEM1 )
      AM1IM = AIMAG( ACOEM1 )
      BM1RE =  REAL( BCOEM1 )
      BM1IM = AIMAG( BCOEM1 )
      T( 4 ) = (2*N - ONE) / ( N*(N - ONE) )
      T( 2 ) = (N - ONE)*(N + ONE) / N
      GQSC   = GQSC + T( 2 )*( AM1RE*ARE + AM1IM*AIM +
     &                         BM1RE*BRE + BM1IM*BIM ) +
     &                T( 4 )*( AM1RE*BM1RE + AM1IM*BM1IM )

      T( 3 ) = 2*N + 1
      QEXT   = QEXT + T( 3 )*( ARE + BRE )
      T( 4 ) = ARE**2 + AIM**2 + BRE**2 + BIM**2
      Qsca   = Qsca + T( 3 )*T( 4 )
      RMM    = - RMM
      SBACK  = SBACK + T( 3 ) * RMM *( ACOE - BCOE )

      T( 2 ) = N *( N + 1 )
      T( 1 ) = T( 3 ) / T( 2 )

      AC = T(1) * ACOE
      BC = T(1) * BCOE
      DO 90 J = 1, NUMANG
         S1( J, 1 ) = S1( J, 1 ) + AC*PI(J, 3) + BC*TAU(J, 3)
         S2( J, 1 ) = S2( J, 1 ) + BC*PI(J, 3) + AC*TAU(J, 3)
   90 CONTINUE

c                               ** Scattering matrix elements for
c                                  supplements of 0-90 degree scattering 
c                                  angles submitted by user
      IF( MOD(N,2).EQ.0 )  THEN

         DO 100 J = 1, NUMANG
            S1( J, 2 ) = S1( J, 2 ) - AC*PI(J, 3) + BC*TAU(J, 3)
            S2( J, 2 ) = S2( J, 2 ) - BC*PI(J, 3) + AC*TAU(J, 3)
  100    CONTINUE

      ELSE

         DO 110 J = 1, NUMANG
            S1( J, 2 ) = S1( J, 2 ) + AC*PI(J, 3) - BC*TAU(J, 3)
            S2( J, 2 ) = S2( J, 2 ) + BC*PI(J, 3) - AC*TAU(J, 3)
  110    CONTINUE

      END IF

c                                      ** Test for convergence of sums
      IF( T(4).GE.1.0E-14 ) THEN

         N  = N + 1
         IF( N.GT.NMX2 )  CALL ERRMSG( 
     &       'MIELAY--Dimensions for W,ACAP not enough. Suggest'//
     &       ' get detailed output, modify routine', .True. )

         DO 120 J = 1, NUMANG

            PI(  J, 1 ) = PI(  J, 2 )
            PI(  J, 2 ) = PI(  J, 3 )
            TAU( J, 1 ) = TAU( J, 2 )
            TAU( J, 2 ) = TAU( J, 3 )

  120    CONTINUE

         ACOEM1  = ACOE
         BCOEM1  = BCOE

         GO TO 70

      END IF

c ***************** End of summing loop ******************


c                            ** Transform complex scattering amplitudes
c                               into elements of real scattering matrix
      DO 140 J = 1, NUMANG

         DO 130 K = 1, 2

            M1( J, K )  = REAL( S1(J,K) )**2 + AIMAG( S1(J,K) )**2
            M2( J, K )  = REAL( S2(J,K) )**2 + AIMAG( S2(J,K) )**2
            S21( J, K ) = REAL(  S1(J,K) )*REAL(  S2(J,K) ) + 
     &                    AIMAG( S1(J,K) )*AIMAG( S2(J,K) )
            D21( J, K ) = AIMAG( S1(J,K) )*REAL( S2(J,K) ) - 
     &                    AIMAG( S2(J,K) )*REAL( S1(J,K) )
  130    CONTINUE

  140 CONTINUE


      T( 1 ) = TWO*RX**2
      QEXT   = T( 1 ) * QEXT
      Qsca   = T( 1 ) * Qsca
      GQSC   = TWO * T( 1 ) * GQSC
      SBACK  = 0.5*SBACK
      QBS    = ( REAL(SBACK)**2 + AIMAG(SBACK)**2 ) / (PINUM*XSHELL**2)

      END

