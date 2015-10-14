      subroutine adapt(f, xknots, coefs, kdimen, ndimen)        
c
c     this algorithm computes a piecewise polynomial approximation
c  of specified smoothness, accuracy and degree.  the input to the
c  computation is
c
c   f      - function being approximated. it must provide values of
c            derivatives up to the order of smoothness specified for
c            the approximation.  the calling sequence is f(x,fderv) and
c            fderv contains the derivatives( see constraint below)
c   a,b    - the endpoints of the interval of approximation
c   accur  - the accuracy required for the approximation
c   smooth - the smoothness required for the approximation
c              = 0  means continuous
c              = 1  means continuous slope
c              = 2  means continuous second derivative, etc.
c   degree - the degree of the polynomial pieces.
c            must have degree gt 2*smooth
c
c     *****  ***** secondary input - items with default values possible
c   charf  - characteristic length of the function f(x). pieces are not
c            longer than this length.
c                   default=(b-a) if degree gt 1, else (b-a)/3
c   norm   - norm to measure the approximation error
c              = 1  l1 approximation (least deviations)
c              = 2  l2 approximation (least squares)
c              = 3  tchebycheff (minimax) approximation
c              =-p  (negative value) general lp approximation
c                   default= 2
c   nbreak - number of special break points in the approximation.
c            associated input variables are
c              xbreak(j)  - location of break points
c              dbreak(j)  - derivative broken at xbreak
c              bleft (j)  - value from left for dbreak derivative
c              bright(j)  -   -    -   right -    -        -
c                   default = 0
c   level  - level of output from adapt
c              = -1 no output at all except for fatal input errors
c              =  0  error conditions noted, final summary
c              =  1  print the approximation out
c              =  2  details of the computation
c              =  3 debug output,  = 4 lots of debug output
c                    default = 0
c   edist  - switch to change from proportional error distribution
c            to fixed distribution. this is primarily of use in
c            approximation of functions with singularities. one should
c            use norm = 1. or so in such cases
c              = 0  proportional distribution
c              = 1  approximate fixed error distribution
c                   attempts to achieve specified accuracy value accur
c              = 2  true fixed error distribution
c
c     the output of the computation consists of 3 parts, each returned
c     to the user in a different way. they are
c
c   xknots,coefs - arrays defining the piecewise polynomial result.
c   coefs    xknots(k)  = knots of the approximation ( k = 1 to knots)
c                         the last one is right end point of interval
c            coefs(k,n) = coefficient of (x - xknot(k))**(n-1) in the
c                         interval xknot(k) to xknot(k+1)
c                           k = 1 to knots-1  and  n = 1 to degree+1
c            these arrays are passed as arguments so as to use variable
c            dimensions.
c               kdimen - dimension used for  xknots in calling program
c               ndimen - coefs is declared  coefs(kdimen,ndimen) in the
c                        calling program.
c               ***** note ***** several small arrays here have fixed
c                        dimensions that limit degree and thus ndimen
c                        should not exceed this limit (currently = 12)
c
c   ppoly  - the piecewise polynomial approximating function.
c            this function subprogram is available to the user at the
c            completion of adapt.
c
c   resulz - a labeled common block with  error,knots  in it
c              knots - number of knots of ppoly
c              error - accuracy of ppoly as estimated by adapt
c
c  ********** dimension constraints **********
c     maxknt - max number of knots taken from user via kdimen
c              arrays with this dimension
c                   coefs   xknots
c     maxpar - max number of parameters per interval ( = 12 currently )
c                   user provided ndimen must have ndimen le maxpar
c                   must have  degree + 1 le maxpar
c              arrays with this dimension (or related values )
c                d       ddtemp  fdervl  fdervr  fdumb   factor
c                fintrp  fleft   fright  powers  xtemp   xintrp  xdd
c              ***** note ***** maxpar also affects argument fderv
c              of function f.  fdervl, fdervr are also involved.
c              should declare fderv of size 6 in f to be safe.
c     maxaux - maximum numbber of auxiliary input( = 20 now ). arrays
c                 xbreak  dbreak  bleft   bright
c     maxstk - max size of active interval stack
c              min interval length is 2**(-maxstk)*(b-a). arrays
c                 xleft   xright
c
c  **********  portability considerations  **********
c
c     this program is in ansi standard fortran.  in addition, it meets
c     all the requirements of the bell labs portable fortran -pfort-
c     except one.  hollerith characters are packed 4/word rather than
c     1/word as specified by pfort.
c     nevertheless, this program is affected in several ways by a
c     change in machine word length and changing to double precision.
c     ***** this version is for the machine with the longest single
c           precision word (cdc).  the length of some constants in
c           the subprogram comput exceeds the capacity of some fortran
c           compilers and can prevent compilation.
c  input-output -- is of three types.
c         fatal error messages - occur in setup,take,put and termin
c                                they cannot be switched off
c         user and debugging output - can be switched off
c     these involve many formats like  e15.8, f12.8, e22.13, etc.
c     some fortran compilers require d-format for double precision.
c     some do not handle e22.13 on machines of shorter word length.
c
c     sumary prints coefficients 5 per line, 6 per line is better
c     for shorter word lengths.  double precision on many machines
c     limits one to 4 per line.
c
c  constants -- the gauss weights and abscissa in compute are given
c               to 15 digits in the comments
c
c  double precision conversion -- requires four steps
c     1. all real variables are declared double precision. this is
c        done by changing real to double precision as all reals are
c        explicitly declared and room is left for this change.  real
c        variables appear before integers in all common blocks.
c
c     2. add d to constants( e.g. 1.0 = 1.0d0 ). adjust length of
c        gauss weights in compute.
c
c     3. change abs,amax1,float  at many places
c
c     4. adjust the interval stack size = dimensions of xleft, xright
c        and value of maxstk.  adjust the value buffer in put to be
c        about .1e-k for a machine with k+2 decimal digits.
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      external f
      real f
c
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
c                      kntdim - kdimen, name changed to put in common
c                      npardm - ndimen, name changed to put in common
      common /resulz/ error, knots
c                      knots = final no. of knots, includes b as one.
c                      error = estimate of error actually achieved.
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
c             kontrl contains generally useful variables
c                      fatal  - switch for detection of fatal error
c                               including excessive interval subdivision
c                               which does not terminate the computation
c                      finish - switch to terminate algorithm
c                      maxs   - see comments earlier
c                      nstack - counter for interval stack, consists of
c                               (xleft(j),xright(j))  j = 1 to nstack
c                      errori - error estimate for top interval
c                      dsctol - tolerance to check discarding intervals
c                      discrd - switch to signal discard of top interval
c                      factor - array of factorials
c                      npar   - number of paremeters = degree + 1
c                      fmesge - string =  **** fatal error *****
c                      interp - number of interior interpolation points
c                               in the normal interval
c                      ibreak - counter on break points
c                      break  - switch for break point in top interval
c                               0     = no break present
c                               left  = break at xleft(nstack)
c                               right = break at xright(nstack)
c                               both  = break at both ends
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c             comdif contains variables used only by comput and friends.
c                      nintrp - number of interior interpolation points
c                               for the current interval
c                      xintrp - interior interpolation points
c                      fintrp - f values at xintrp points
c                      leftx  - multiplicity of interpolation at xleft
c                               = no. of derivatives matched at xleft
c                      fleft  - values of f and its derivatives at xleft
c                      rightx - multiplicity of interpolation at xright
c                      fright - values of f and derivatives at xright
c                      ddtemp - the array of divided differences
c                      xdd    - the x values for ddtemp with proper
c                               multiplicities of xleft and xright
      common /saveit/ iknot
c
c------------------------ main control program -------------------------
c
c   ***** note - arguments below are for readability only    *****
c   *****        except for f and xknots,coefs,kdimen,ndimen *****
c
c                   check input, initial computations, print problem
c
      call setup(xknots, coefs, kdimen, ndimen)
c
c         check for fatal error in problem specification
      if (fatal) return
c
c                   loop over processing of intervals
   10 call take(interv)
c
      call comput(f, approx, interv)
c
c                   check for discarding intervals
      call check(funct, charct)
c
c            put new intervals on stack or discard, update status
      call put(interv, xknots, coefs, kdimen, ndimen)
c
      call termin(test, and, print, xknots, kdimen)
c
      if (.not.finish) go to 10
c
      call sumary(xknots, coefs, kdimen, ndimen)
      return
      end
      subroutine check(funct, char)   
c
c        ===============================================================
c
c **  this program checks for discarding interval, applies various
c     tests about discarding involving edist and charf.
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real dtest, dx
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c
      discrd = .false.
      dx = xright(nstack) - xleft(nstack)
c
c         compute dtest for implementing various types of adaptive apprx
      if (edist.eq.0) dtest = dx*dsctol
c         for the approximate fixed error distribution type we estimate
c         the final number of knots by( limiting it a little )
c             (nstack+knots+2)((b-a)/(xright-a))
      if (edist.eq.1) dtest = dsctol/(float(nstack+knots+2)*
     *     min((b-a)/(xright(nstack)-a),5.d0))
      if (edist.eq.2 .or. norm.eq.3.d0) dtest = dsctol
c
c             check for discard of interval
      if (errori.le.dtest) discrd = .true.
c
c         check characteristic length of function
      if (dx.ge.charf) discrd = .false.
      return
      end
      subroutine comput(f, approx, interv)  
c
c        ===============================================================
c
c **  this program computes the piecewise polynomial approximation on
c     the current interval. it also estimates the error
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real absc, dx, errint, f, fdervl, fdervr, fdumb, polydd, wgts
      dimension absc(4), wgts(4), fdervl(5), fdervr(5), fdumb(5)
      external f, polydd
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
      equivalence (fleft(2),fdervl(1)), (fright(2),fdervr(1))
c     fifteen digit values for these gauss integration constants are
c  .861136311594053 .339981043584856 .347854845137454 .652145154862546
      data absc(1) /-.86113631159405d0/
      data absc(2) /-.33998104358486d0/
      data absc(3) /.33998104358486d0/
      data absc(4) /.86113631159405d0/
      data wgts(1) /.34785484513745d0/
      data wgts(2) /.65214515486255d0/
      data wgts(3) /.65214515486255d0/
      data wgts(4) /.34785484513745d0/
c
c             compute interpolation information
      nintrp = degree - 2*smooth - 1
c
c          increase number of interpolation points if break points are
c          specified with fewer derivatives than smooth
      if (break.eq.left .or. break.eq.right) nintrp = nintrp + smooth -
     * dbreak(ibreak)
      if (break.eq.both) nintrp = nintrp + 2*smooth - dbreak(ibreak) -
     * dbreak(ibreak+1)
      if (nintrp.eq.0) go to 20
c             generate equal spaced interpolation points
      dx = (xright(nstack)-xleft(nstack))/float(nintrp+1)
      do 10 j=1,nintrp
        xintrp(j) = xleft(nstack) + float(j)*dx
   10 continue
c
c             get left and right f-values, put f-value in first element
c             of arrays fleft and fright.  get derivatives back as
c             other elements via the subarrays fdervl and fdervr.
   20 fleft(1) = f(xleft(nstack),fdervl)
      fright(1) = f(xright(nstack),fdervr)
      leftx = smooth + 1
      rightx = leftx
c            get f-values at other interpolation points, if any
      if (nintrp.eq.0) go to 40
      do 30 j=1,nintrp
        fintrp(j) = f(xintrp(j),fdumb)
   30 continue
c
c          check for break points, modify values if necessary
   40 continue
      if (break.ne.left) go to 50
      leftx = dbreak(ibreak) + 1
      fleft(leftx) = bright(ibreak)
   50 if (break.ne.right) go to 60
      rightx = dbreak(ibreak) + 1
      fright(rightx) = bleft(ibreak)
   60 if (break.ne.both) go to 70
      leftx = dbreak(ibreak) + 1
      rightx = dbreak(ibreak+1) + 1
      fleft(leftx) = bright(ibreak)
      fright(rightx) = bleft(ibreak+1)
   70 continue
c
c           compute divided differences, newton form of polynomial
      call newton(leftx, rightx, nintrp)
c
c         compute norm of error of this appromimation using four pts
c         add 50 percent fudge factor
      errori = errint(f,polydd,xleft(nstack),xright(nstack),absc,wgts)
      errori = 1.5d0*errori
      return
      end
      real function errint(f, fit, aaa, bbb, points, weight)       
c
c        ===============================================================
c
c **  this function does a four point integration rule for the
c     absolute value of the difference of two functions( f and fit )
c                   abs( f(x) - fit(x) )**norm
c     the integration uses the points and weights given and scaled
c     from (-1,1) to (aaa,bbb)
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real aaa, abmid, ba, bbb, f, fdumb, fit, p, pj, points, weight
      real er, f1, f2
      dimension fdumb(5), points(4), weight(4)
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c         compute midpoint = abmid and half length = ba of interval
      abmid = (aaa+bbb)/2.d0
      ba = (bbb-aaa)/2.d0
      pj = abmid + ba*points(1)
c
c         test for tchebycheff (minimax) norm which uses special code
      if (norm.eq.3.d0) go to 20
c
c         have general lp norm or least squares or least deviations
      p = abs(norm)
c              initialize the quadrature rule
      errint = abs(f(pj,fdumb)-fit(pj))**p*weight(1)
c              loop through remaining points
      do 10 j=2,4
        pj = abmid + ba*points(j)
c             debug debug debug start of the computation
        f1 = f(pj,fdumb)
        f2 = fit(pj)
        er = abs(f1-f2)**p
        if (level.ge.4 .and. (knots.eq.1 .or. finish)) write (6,99999)
     *   pj, f1, f2, er
        errint = errint + abs(f(pj,fdumb)-fit(pj))**p*weight(j)
   10 continue
      errint = errint*ba
      go to 40
c
c         tchebycheff norm
   20 continue
c             find max error on points
c               initialize
      errint = abs(f(pj,fdumb)-fit(pj))
c              loop through the remaining points
      do 30 j=2,4
        pj = abmid + ba*points(j)
        errint = max1(errint,abs(f(pj,fdumb)-fit(pj)))
   30 continue
   40 continue
c             debug  debug  debug  debug
      if (level.ge.3) write (6,99998) errint, aaa, bbb
      return
99999 format (5(3h --), 31hdebug error curve,pj,f1,f2,er =, 4f15.8)
99998 format (15x, 9herrint = , f20.15, 4h on , 2f15.8)
      end
      subroutine newton(nl, nr, ni)     
c
c        ===============================================================
c
c **  this program computes the divided differences array as follows
c         nl coalesced points on left   - deriv values in fleft
c         nr coalesced points on right  -   -      -    - fright
c         ni distinct  points inbetween - fnctn    -    - fintrp
c
c     the points are ordered xl = xleft (nstack)
c                            xr = xright(nstack)
c                            xintrp array
c
c         layout of the ddtemp divided difference array
c
c     nl=6    llllll****ii
c     nr=4    lllll****ii     l = first triangle
c     ni=2    llll****ii
c             lll****ii       r = second triangle
c             ll****ii
c             l****ii         * = fill between triangles
c             rrrrii
c             rrrii           i = completion for interpolation points
c             rrii
c             rii       idif = horizontal coord. = difference order
c             ii        ipt  = vertical coord. associated with points
c             i
c
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real difff, diffx
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c                 main calculation of divided differences
c         define a few short constants
      nl1 = nl - 1
      nl2 = nl + 1
      nr1 = nr - 1
      nr2 = nr + 1
      nrl = nr + nl
c
c         put x-values in a single array with nddx = nl+nr+ni points
      do 10 nddx=1,nl
        xdd(nddx) = xleft(nstack)
   10 continue
      nddx = nl
      do 20 k=1,nr
        nddx = nddx + 1
        xdd(nddx) = xright(nstack)
   20 continue
c            check if there are any interpolation points to add to xdd
      if (ni.eq.0) go to 40
      do 30 k=1,ni
        nddx = nddx + 1
        xdd(nddx) = xintrp(k)
   30 continue
c
c           fill border of first triangle - size nl.
   40 continue
c         top border
      do 50 idif=1,nl
        ddtemp(idif,1) = fleft(idif)/factor(idif)
   50 continue
      if (nl1.eq.0) go to 70
c                   bottom border
      do 60 idif=1,nl1
        ipt = nl2 - idif
        ddtemp(idif,ipt) = ddtemp(idif,1)
   60 continue
c
c          fill border of second triangle - size nr
   70 continue
c         top border
      do 80 idif=1,nr
        ddtemp(idif,nl2) = fright(idif)/factor(idif)
   80 continue
      if (nrl.eq.nl2) go to 100
c                   bottom border
      do 90 idif=1,nr1
        ipt = nrl + 1 - idif
        ddtemp(idif,ipt) = ddtemp(idif,nl2)
   90 continue
c
c           fill parallogram between the two triangles just filled
c          fill entries parallel to bottom of first triangle
  100 continue
c
c         loop stepping along top side of second triangle
      do 120 j=2,nr2
        idif = j
c             loop stepping parallel to bottom side of first triangle
        do 110 k=2,nl2
          ipt = nl + 2 - k
          difff = ddtemp(idif-1,ipt+1) - ddtemp(idif-1,ipt)
          ipt2 = ipt - 1 + idif
          diffx = xdd(ipt2) - xdd(ipt)
          ddtemp(idif,ipt) = difff/diffx
          idif = idif + 1
  110   continue
  120 continue
c             debug  debug  debug  debug
      if (level.ge.4 .and. knots.le.1) write (6,99999) nr2, nl2, idif,
     * ipt, difff, diffx
c
c         fill in bottom diagonals for interpolation points, if any
      if (ni.eq.0) go to 150
c         loop through the interpolatation points
      do 140 j=1,ni
        idif = 2
        nrlj = nrl + j
        ddtemp(1,nrlj) = fintrp(j)
c         loop through the differences (idif index)
        nrlj1 = nrlj - 1
        do 130 k=1,nrlj1
          ipt = nrlj - k
          difff = ddtemp(idif-1,ipt+1) - ddtemp(idif-1,ipt)
          diffx = xdd(nrlj) - xdd(ipt)
          ddtemp(idif,ipt) = difff/diffx
          idif = idif + 1
  130   continue
  140 continue
  150 continue
      return
99999 format (9x, 30hnr2,nl2,idif,ipt,difff,diffx =, 4i3, 2f12.6)
      end
      real function polydd(x)                                
c
c        ===============================================================
c
c **  this function evaluates the current polynomial piece represented
c     by the divided differences ddtemp on the point set xdd.
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real x
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c
      polydd = ddtemp(degree+1,1)
      do 10 k=1,degree
        j = degree + 1 - k
        polydd = ddtemp(j,1) + (x-xdd(j))*polydd
   10 continue
      return
      end
      subroutine ppfit4(f, xlft, xrgt, epsln, npiece, errest, xknots,  
     * coefs, kdimen, ndimen, ndeg, nsmth, emeas, lprnt, foscil, atype,
     * kbreak, brakpt, kdervb, vallft, valrgt)
c
c        ===============================================================
c
c **  this set of four control programs set varying numbers of default
c     values for arguments.  it uses entry statements which are done
c     differently by various fortrans.  for this reason entry statements
c     are only indicated by comment cards.  writing four separate rou-
c     tines approximately triples the length of the total code.
c
c     the following tabulates the internal and external names of the
c     arguments along with their default values for the various ppfit.
c     note that this allows the arguments to be put into a common block
c     and avoids long argument lists within adapt itself.
c
c
c    internal   external     default value    set by ppfit number
c     f          f              none
c     a,b        a,b            none
c     accur      epsln          none
c     knots      npiece         output
c     error      errest         output
c     xknots     xknots         output
c     coefs      coefs          output
c     kdimen     kdimen         none
c     ndimen     ndimen         none
c     degree     ndeg           3                1
c     smooth     nsmth          0                1
c     norm       emeas          3                1 2
c     level      lprnt          1                1 2
c     charf      foscil         variable         1 2
c     edist      atype          1                1 2 3
c     nbreak     kbreak         0                1 2 3
c     xbreak     brakpt         -
c     dbreak     kdervb         -
c     bleft      vallft         -
c     bright     valrgt         -
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      real brakpt, emeas, epsln, errest, foscil, f, vallft, valrgt,
     * xlft, xrgt
      dimension brakpt(20), kdervb(20), vallft(20), valrgt(20)
      integer atype
      external f
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
      common /saveit/ iknot
c
c
c             skip all default setting for ppfit4
      go to 10
c
c                 set defaults for ppfit1
c
c         ****** since entry is non-standard  ******
c         ****** the natural way to implement ******
c         ****** the other ppfits is only     ******
c         ****** indicated in the comments    ******
c
c     entry ppfit1
c     entry ppfit1(f,xlft,xrgt,epsln,npiece,errest,xknots,coefs,
c    a             kdimen,ndimen)
      ndeg = 3
      nsmth = 0
c
c                 set defaults for ppfit2
c     entry ppfit2
c     entry ppfit2(f,xlft,xrgt,epsln,npiece,errest,xknots,coefs,
c    a             kdimen,ndimen,ndeg,nsmth)
      emeas = 3.0d0
      lprnt = 1
      foscil = b - a
      if (ndeg.eq.2) foscil = .5d0*foscil
      if (ndeg.le.1) foscil = (b-a)/3.d0
c
c                 set defaults for ppfit3
c     entry ppfit3
c     entry ppfit3(f,xlft,xrgt,epsln,npiece,errest,xknots,coefs,
c    a             kdimen,ndimen,ndeg,nsmth,emeas,lprnt,foscil)
      atype = 1
      ksing = 0
      kbreak = 0
c
c         put input into common inputz by changing to internal names
   10 a = xlft
      b = xrgt
      accur = epsln
      degree = ndeg
      smooth = nsmth
      norm = emeas
      level = lprnt
      charf = foscil
      edist = atype
      nbreak = kbreak
      if (nbreak.le.0 .or. nbreak.ge.21) go to 30
      do 20 k=1,nbreak
        xbreak(k) = brakpt(k)
        dbreak(k) = kdervb(k)
        bleft(k) = vallft(k)
        bright(k) = valrgt(k)
   20 continue
   30 continue
c
      call adapt(f, xknots, coefs, kdimen, ndimen)
      npiece = knots
      errest = error
      return
      end
      real function ppoly(t, xknots, coefs, kdimen, ndimen)   
c
c        ===============================================================
c
c **  this function evaluates the piecewise polynomial approximation
c     computed in adapt
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      real t, x
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
      common /saveit/ iknot
c           start search for right interval from point of last
c           evaluation of ppoly
c
c           forward searching loop
   10 if (t.lt.xknots(iknot+1)) go to 20
      iknot = iknot + 1
      if (iknot.lt.knots) go to 10
      iknot = knots - 1
c
c                reached right end of interval
      go to 30
c
c           backward searching loop
   20 if (t.ge.xknots(iknot)) go to 30
      iknot = iknot - 1
      if (iknot.gt.1) go to 20
      iknot = 1
c
c                reached left end of interval
c
c                evaluate from powers based at xknot(iknot)
c                use nested multiplication
   30 x = t - xknots(iknot)
      ppoly = coefs(iknot,npar)
      do 40 k=1,degree
        kk = npar - k
        ppoly = coefs(iknot,kk) + x*ppoly
   40 continue
      return
      end
      subroutine ptrans(d, powers)
c
c        ===============================================================
c
c **  this program converts polynomial representation from divided
c     difference to power form.  there are coalesced points on each
c     end of the interval (xl,xr) = (xleft(nstack),xright(nstack)).
c     the number coalesced at each end is leftx and rightx.
c     and there are nintrp other pts xintrp(k)  inbetween them.
c     see subroutine newton for more details
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real d, powers, shift, xl, xr, xtemp
      dimension d(12,12), powers(12), xtemp(12)
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c         set some short local variable names
      xl = xleft(nstack)
      xr = xright(nstack)
      nl = leftx
      nl1 = nl + 1
      nr = rightx
      ni = nintrp
      nrl = nr + nl
      nri = nr + ni
      nri1 = nri - 1
      nrli = nrl + ni
c
c         starting representation is (assuming xl = 0 )
c
c    d(1) +d(2)x +d(3)x**2 + --- +d(nl)x**(nl-1)
c   +(x**nl)*( d(nl+1)(+d(nl+2)(x-xr)**2 + --- +d(nl+nr)*(x-xr)**(nr-1)
c        *((x-xr)**nr)*(d(nl+nr+1) + d(nl+nr+2)*(x-xintrp(1))
c                      +d(nl+nr+3)*(x-xintrp(1))(x-xintrp(2)) + ---))
c
c         strategy is to first convert the part from the interp. pts.
c         to poly in (x-xr).  this poly then has origin shifted to xl.
c
c     the conversion of the interp part is done explicitly for degree
c     two or less and done by synthetic division for higher degrees
c
c   d1 + d2(x-x1) +d3(x**2-(x1+x2)x +x1*x2)
c
c     the resulting coefficients are put in the array powers
c
c             debug  debug  debug  debug
      if (level.eq.4) write (6,99999) (d(k,1),k=1,npar)
      if (ni.eq.0) go to 100
c             build up the polynomial for the interpolation points
c
c         use special formulas for ni less than 3
      if (ni.eq.1) go to 10
      if (ni.eq.2) go to 20
      go to 30
   10 powers(1) = d(nrl+1,1)
      go to 80
   20 powers(1) = d(nrl+1,1) + (xr-xintrp(1))*d(nrl+2,1)
      powers(2) = d(nrl+2,1)
      go to 80
c
c         conversion by repeated synthetic division
   30 ni1 = ni - 1
c         initialize the powers and xtemp arrays
      do 40 k=1,ni
        xtemp(k) = xintrp(k)
        nrlk = nrl + k
        powers(k) = d(nrlk,1)
   40 continue
c
c         do the repeated synthetic division to replace the xtemp
c         = xintrp points of the newton expansion by the xr points
      do 70 k=1,ni1
c              powers(ni) is fixed and set above
        do 50 ii=1,ni1
          i = ni - ii
          powers(i) = powers(i) + (xr-xtemp(i))*powers(i+1)
   50   continue
c              shift the newton expansion pts. up, put in one more xr
        do 60 ii=1,ni1
          i = ni - ii
          xtemp(i+1) = xtemp(i)
   60   continue
        xtemp(1) = xr
   70 continue
      if (level.eq.4) write (6,99998) (powers(k),k=1,ni)
   80 continue
c             shift the coefficients to the top of the powers array
      do 90 k=1,ni
        l = ni + 1 - k
        ltop = l + nrl
        powers(ltop) = powers(l)
   90 continue
c
c             have the interpolation pt. coefs. in the array powers
  100 continue
c             put the remaining divided diffs into the powers array
      do 110 j=1,nrl
        powers(j) = d(j,1)
  110 continue
c
c             debug  debug  debug  debug
      if (level.eq.4) write (6,99997) (powers(k),k=nl1,npar)
c         transform the origin of the polynomial from xr to xl
c         we use repeated synthetic division
      if (nri.eq.1) go to 140
      shift = xr - xl
      khi = nri1
c             loop through the coefficients
      do 130 j=2,nri
c                   synthetic division loop
        do 120 k=1,khi
          koef = nrli - k
          powers(koef) = powers(koef) - shift*powers(koef+1)
  120   continue
        khi = khi - 1
  130 continue
  140 continue
c     the coefficients are now of the power form with origin xl
c             debug  debug  debug  debug
      if (level.eq.4) write (6,99996) (powers(k),k=1,npar)
      return
99999 format (15x, 26hdebug ptrans, orig input =/(5x, 8e15.5))
99998 format (10x, 17hinterp part coefs/(5x, 8e15.5))
99997 format (10x, 25hright + interp part coefs/(5x, 8e15.5))
99996 format (10x, 18hfinal coefficients/(5x, 8e15.5))
      end
      subroutine put(interv, xknots, coefs, kdimen, ndimen)       
c
c        ===============================================================
c **  this program puts intervals on the stack or discards them.
c     when an interval is discarded a new knot is found. then this
c     program updates the error estimate, the xknot array, transforms
c     the polynomial to the power form and put the coefficients into
c     the array coefs.  it also checks for passing break points
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      real buffer, dx, powers, p, ratio
      dimension powers(12)
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
      data buffer /1.d-12/
c            check for discarding the interval
      if (discrd) go to 30
c         subdivide interval and place on stack
      if (nstack.lt.maxstk) go to 10
c             fatal error, exceeded active stack size
      fatal = .true.
      if (level.ge.0) write (6,99999) fmesge, maxstk, xleft(nstack),
     * xright(nstack)
      discrd = .true.
      go to 30
c
   10 dx = (xright(nstack)-xleft(nstack))*.5d0
c              check for small intervals
      ratio = dx/(abs(a)+abs(b))
      if (ratio.gt.buffer) go to 20
      discrd = .true.
      fatal = .true.
      if (level.ge.0) write (6,99998) xleft(nstack), xright(nstack)
      go to 30
   20 continue
      nstack = nstack + 1
      xleft(nstack) = xleft(nstack-1)
      xleft(nstack-1) = xright(nstack-1) - dx
      xright(nstack) = xleft(nstack-1)
c             debug  debug  debug  debug
      if (level.ge.3) write (6,99997) nstack, xleft(nstack),
     * xright(nstack)
      return
c
c            discard interval, update global error, xknots and coefs
   30 continue
c
      p = abs(norm)
      if (norm.eq.3.d0) error = max1(error,errori)
      if (norm.ne.3.d0) error = (error**p+errori)**(1.d0/p)
c
c              check for passing break points
      if (break.eq.left .or. break.eq.both) ibreak = ibreak + 1
c             debug  debug  debug  debug
      if (level.ge.3) write (6,99996) nstack
c
c             transform representation of polynomial from divided
c             differences to powers of x with origin at xknots (knots)
      call ptrans(ddtemp, powers)
c
c             put coefs into the main array
      do 40 k=1,npar
        coefs(knots,k) = powers(k)
   40 continue
c            put the new knots in xknots
      knots = knots + 1
      xknots(knots) = xright(nstack)
      nstack = nstack - 1
      return
99999 format (//2x, 6a4, 36hinterval divided too much, exceeded ,
     * 6hlimit , i3, 22h on interval stack at /20x, 2e18.8/2x, 6hinterv,
     * 38hal discarded and computation continued)
99998 format (25x, 23hgot short interval ****, 2e22.13, 12h **** discar,
     * 4hd it)
99997 format (15x, 13hput interval , i3, 10h on stack , 2f12.8)
99996 format (15x, 16hdiscard interval, i5)
      end
      subroutine setup(xknots, coefs, kdimen, ndimen)           
c
c        ===============================================================
c
c **  this program checks the input data and initializes the computation
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      integer hmsge(6), nmsge(4,4)
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c              iknot is a counter used in the output approx. ppoly
      common /saveit/ iknot
      data hmsge(1), hmsge(2), hmsge(3), hmsge(4), hmsge(5), hmsge(6) /
     * 4h ***,4h* fa,4htal ,4herro,4hr **,4h**  /
      data nmsge(1,1), nmsge(1,2), nmsge(1,3), nmsge(1,4) /4hleas,
     * 4ht de,4hviat,4hions/
      data nmsge(2,1), nmsge(2,2), nmsge(2,3), nmsge(2,4) /4hleas,
     * 4ht sq,4huare,4hs   /
      data nmsge(3,1), nmsge(3,2), nmsge(3,3), nmsge(3,4) /4hmini,
     * 4hmax ,4hnorm,4h    /
      data nmsge(4,1), nmsge(4,2), nmsge(4,3), nmsge(4,4) /4hgene,
     * 4hral ,4hl-p ,4hnorm/
      data kleft, kright, kboth /4hleft,4hrite,4hboth/
c
c         put data statement items into common variables
c         various compilers allow more efficient ways
c
      left = kleft
      right = kright
      both = kboth
      do 10 k=1,6
        fmesge(k) = hmsge(k)
   10 continue
c -------- set current values of limits on dimensions ------------------
      kntdim = kdimen
      npardm = ndimen
      maxknt = kntdim
      maxstk = 50
      maxpar = min0(12,npardm)
      maxaux = 20
c -------- check input data --------------------------------------------
      fatal = .false.
      if ((level+1)*level*(level-1)*(level-2).eq.0) go to 20
c         fixed error on print control level
      write (6,99999) level
c         for debugging do not change output level
c         for production version set level = 0
   20 if (a.lt.b) go to 30
c         fatal error in interval (a,b)
      fatal = .true.
      write (6,99998) fmesge, a, b
   30 if (accur.gt.0.0d0) go to 40
c         fatal error, negative accuracy
      fatal = .true.
      write (6,99997) fmesge, accur
   40 if (degree.lt.maxpar) go to 50
c         fatal error, degree exceeds maximum allowed value
      if (level.ge.0) write (6,99996) fmesge, degree
      fatal = .true.
   50 if (2*smooth.lt.degree) go to 60
c         fatal error, smooth and degree incompatible
      fatal = .true.
      write (6,99995) fmesge, smooth, degree
   60 if (charf.gt.(b-a)/float(maxknt)) go to 70
c         fatal error, charf is too small
      fatal = .true.
      write (6,99994) fmesge, charf
   70 if ((norm-1.d0)*(norm-2.d0)*(norm-3.d0).eq.0.0d0
     $     .or. norm.lt.0.0d0) go to 80
c         fatal error, undefined norm requested
      if (level.ge.0) write (6,99993) fmesge, norm
      fatal = .true.
   80 if (nbreak.ge.0 .and. nbreak.le.maxaux) go to 90
c         fatal error in number of break points
      fatal = .true.
      write (6,99992) fmesge, nbreak
c
   90 if (nbreak.eq.0) go to 150
c         check the break point data, monotonicity and degree
      j = 1
      if (xbreak(1).lt.a .or. xbreak(nbreak).gt.b) go to 110
      if (nbreak.eq.1) go to 120
      do 100 j=2,nbreak
        if (xbreak(j-1).ge.xbreak(j)) go to 110
  100 continue
      go to 120
  110 fatal = .true.
c         break points are not monotonic
      write (6,99991) fmesge, xbreak(j)
  120 limsm = (degree-1)/2
      do 130 j=1,nbreak
        if (dbreak(j).lt.0 .or. dbreak(j).gt.limsm) go to 140
  130 continue
      go to 150
c         bad value in derivative breaks
  140 fatal = .true.
      write (6,99990) fmesge, dbreak(j)
c
  150 if (edist*(edist-1)*(edist-2).eq.0) go to 160
c         fatal error, illegal error distribution requested
      if (level.ge.0) write (6,99989) fmesge, edist
      fatal = .true.
  160 continue
c
c -------- initialization of variables ---------------------------------
c
c             active interval stack
      nstack = 1
      xleft(1) = a
      xright(1) = b
c             termination and error values
      finish = .false.
      error = 0.d0
      dsctol = accur**abs(norm)
      if (edist.eq.0) dsctol = dsctol/(b-a)
      if (norm.eq.3.d0) dsctol = accur
c             miscellaneous variables and pointers
      ibreak = 1
      knots = 1
      iknot = 1
      interp = degree + 2 - 2*smooth
      xknots(1) = a
      npar = degree + 1
c
c             compute array of npar factorials, factor(k)= k-1 factorial
      factor(1) = 1.d0
      factor(2) = 1.d0
      do 170 k=3,npar
        factor(k) = float(k-1)*factor(k-1)
  170 continue
c
c -------- print out of problem to be solved ---------------------------
      if (level.le.0) go to 180
      nmes = norm
      if (norm.lt.0.0d0) nmes = 4
      write (6,99988) a, b, degree, smooth, accur, (nmsge(nmes,j),j=1,4)
      write (6,99987) charf, norm
      if (edist.eq.0) write (6,99986)
      if (edist.eq.1) write (6,99985)
      if (edist.eq.2) write (6,99984)
      if (edist.eq.-1) write (6,99983)
      if (nbreak.eq.0) go to 180
      write (6,99982) nbreak
      write (6,99981) (xbreak(j),dbreak(j),bleft(j),bright(j),j=1,
     * nbreak)
  180 continue
      return
99999 format (//2x, 29hillegal output level control , i2, 9h set to 0)
99998 format (//2x, 6a4, 20h incorrect interval , 2e15.5)
99997 format (//2x, 6a4, 28h illegal accuracy requested , e15.5)
99996 format (//2x, 6a4, 6hdegree, i3, 29h exceeds maximum allowed valu,
     * 1he)
99995 format (//2x, 6a4, 35h incompatible degree and smoothness, 2i5)
99994 format (//2x, 6a4, 40h characteristic oscillation length of f ,
     * e15.5, 14h is too small )
99993 format (//2x, 6a4, 20h error measure norm , e15.5, 11h not allowe,
     * 1hd)
99992 format (//2x, 6a4, 11h in number , i4, 17h of break points )
99991 format (//2x, 6a4, 30h break points not in order at , e15.5)
99990 format (//2x, 6a4, 20h illegal derivative , i3, 10h at break )
99989 format (//2x, 6a4, 32hillegal error distribution type , i3)
99988 format (//2x, 45h++++++ piecewise polynomial approximation on ,
     * 9hinterval , 2e15.4, 11h of degree , i2, 6h with , i2, 7h contin,
     * 17huous derivatives /10x, 23h accuracy requested is , e15.4,
     * 13h measured by , 4a4)
99987 format (10x, 43hother input/default variables are foscil = ,
     * e15.4, 10h, emeas = , e15.4)
99986 format (15x, 9(2h--), 33h proportional error distribution )
99985 format (15x, 9(2h--), 36happroximate fixed error distribution)
99984 format (15x, 9(2h--), 25h fixed error distribution)
99983 format (15x, 9(2h--), 28hmodified proportional error , 8hdistribu,
     * 4htion)
99982 format (i12, 36h break points specified with data = , 9h location,
     * 32h, derivative, left+right values )
99981 format (5x, 2(e20.5, i4, e16.5, e15.5))
      end
      subroutine sumary(xknots, coefs, kdimen, ndimen)           
c
c        ===============================================================
c
c **  this program prints out a summary of results of adapt
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
      if (level.eq.-1) return
c                   level 0 output
      write (6,99999) degree, smooth, knots, error
      if (level.eq.0) return
c
c                   level 1 output - knots and coefficients
      knots1 = knots - 1
      if (knots.gt.15) go to 20
c
c               format for small no. of knots
      write (6,99998)
      do 10 k=1,knots1
        write (6,99997) k, xknots(k), coefs(k,1)
        write (6,99996) (coefs(k,j),j=2,npar)
   10 continue
      go to 40
c
c               format for lots of knots
   20 continue
      write (6,99995)
      do 30 k=1,knots1
        write (6,99994) k, xknots(k), (coefs(k,j),j=1,npar)
   30 continue
   40 return
99999 format (///48h --- adaptive piecewise polynomial approximation,
     * 11h of degree , i2, 6h with , i2, 12h continuous , 10hderivative,
     * 8hs needed, i4, 17h knots for error , e10.4)
99998 format (8x, 13hknot location, 13x, 21hx-power coefficients ,
     * 27hrelative to knot locations )
99997 format (i10, 2e20.12)
99996 format (30x, e20.12)
99995 format (3x, 39hk   k-th interior knot     powers of x , 7hrelativ,
     * 20he to knot locations )
99994 format (i5, e25.12, 4e22.12/(e30.12, 4e22.12))
      end
      subroutine take(interv)                            
c
c        ===============================================================
c
c **  this program takes an active interval off the top of the stack
c     it also does most of the work of locating and handling
c     break points
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real dx, ratio, buffer
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
      data buffer /1.d-12/
c
c         check for break point
      break = 0
      if (nbreak.eq.0 .or. ibreak.gt.nbreak) go to 20
      if (xbreak(ibreak).gt.xright(nstack)) go to 20
c
c               set control variable break, check for location
      if (xbreak(ibreak).gt.xleft(nstack)) go to 10
      break = left
      if (ibreak.eq.nbreak) go to 20
c             check for second break point in this interval
      if (xbreak(ibreak+1).ge.xright(nstack)) go to 20
c                 next break is inside interval, split top interval
      break = both
c                  check exceeding stack limit. if so, stop
      if (nstack.eq.maxstk) go to 30
c                     dont split very small intervals, stop instead
      dx = xbreak(ibreak+1) - xleft(nstack)
      ratio = dx/(abs(a)+abs(b))
      if (ratio.le.buffer) go to 40
      nstack = nstack + 1
      xleft(nstack) = xleft(nstack-1)
      xright(nstack) = xbreak(ibreak+1)
      xleft(nstack-1) = xright(nstack)
      if (level.ge.2) write (6,99999) nstack, xleft(nstack),
     * xright(nstack)
      go to 20
c
   10 break = right
c                 check to see if break is already at right end point
      if (xbreak(ibreak).ge.xright(nstack)) go to 20
c                 the break is inside interval, split top interval
c                  check exceeding stack limit. if so, stop
      if (nstack.eq.maxstk) go to 30
c                     dont split very small intervals, stop instead
      dx = xbreak(ibreak) - xleft(nstack)
      ratio = dx/(abs(a)+abs(b))
      if (ratio.le.buffer) go to 40
      nstack = nstack + 1
      xleft(nstack) = xleft(nstack-1)
      xright(nstack) = xbreak(ibreak)
      xleft(nstack-1) = xright(nstack)
      if (level.ge.2) write (6,99998) nstack, xleft(nstack),
     * xright(nstack)
   20 continue
c             debug  debug  debug
      if (level.ge.3) write (6,99997) xleft(nstack), xright(nstack),
     * break
      return
c
c        unable to proceed because the stack limit is reached with
c        multiple breakpoints inside the smallest one.
   30 fatal = .true.
      finish = .true.
      if (level.ge.0) write (6,99996) fmesge, maxstk, xleft(nstack),
     * xright(nstack), xbreak(ibreak)
      return
c
c         splitting intervals to accomodate break points has led to
c         an excessively small interval
   40 fatal = .true.
      finish = .true.
      if (level.ge.0) write (6,99995) fmesge, xleft(nstack),
     * xbreak(ibreak), xright(nstack)
      return
99999 format (10x, 28h++++ split top interval, get, i3, 2e18.8,
     * 17h for break = both)
99998 format (10x, 28h++++ split top interval, get, i3, 2e18.8,
     * 19h for break = right )
99997 format (15x, 14htake interval , 2f12.8, 11h    break =, a4)
99996 format (//2x, 6a4, 42h interval divided too much, exceeded limit,
     * i3, 21h on interval stack at/20x, 2e18.8/2x, 12hbreak point ,
     * e18.8, 48h present requires splitting interval, stop adapt)
99995 format (//2x, 6a4, 43hbreak point adjustment has led to a fatally,
     * 1h , 40hsmall interval, the points involved are /20x, 3e18.8)
      end
      subroutine termin(test, and, print, xknots, kdimen)        
c
c        ===============================================================
c
c **  this program tests for termination and gives intermediate output
c
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     * errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     * xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     * xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     * rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen)
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     * dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
      common /resulz/ error, knots
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
c
c         intermediate output - according to levels
c         no intermediate output for levels -1,0,1
      if (level.le.1) go to 40
c
c             level 2 output - only for discarded interval
      if (.not.discrd) go to 10
      write (6,99999) knots, xknots(knots), errori, error
      if (break.eq.left) write (6,99998)
      if (break.eq.right) write (6,99997)
      if (break.eq.both) write (6,99996)
c
      go to 20
c             debug output  -  level 3
   10 if (level.eq.2) go to 40
c
c                  interval summary
      write (6,99995) nstack, xleft(nstack), xright(nstack), errori
      if (break.ne.0) write (6,99994) ibreak, xbreak(ibreak)
c
c             debug polynomial operations  -  level 4
   20 continue
      if (level.le.3) go to 40
      write (6,99993) leftx, rightx, nintrp
      do 30 k=1,npar
        knpar = npar + 1 - k
        write (6,99992) k, (ddtemp(i,k),i=1,knpar)
   30 continue
   40 continue
c             test for normal termination
      if (nstack.eq.0) go to 50
c
c             test for abnormal termination
      if (knots.lt.maxknt) return
c
c            have exceeded limit on knots
      write (6,99991) fmesge, maxknt, xknots(knots), error
      fatal = .true.
c
c         terminate computation
   50 finish = .true.
      return
99999 format (2x, 9h**** knot, i4, 4h at , e16.5, 17h, with local and ,
     * 16hglobal errors = , 2e15.4)
99998 format (20x, 32hbreak point on left of interval )
99997 format (20x, 33hbreak point on right of interval )
99996 format (20x, 37hbreak point at both ends of interval )
99995 format (10x, 19h++++ stack element , i3, 3h = , e20.5, e15.5,
     * 19h, with local error , e15.4, 17h placed on stack )
99994 format (15x, 14hbbbreak point , i3, 3h = , e15.4, 9h involved)
99993 format (5x, 44h---- divided difference info ---- nl,nr,ni =, 3i3)
99992 format (5x, 12hdiv diff row, i3, 9f12.5/93x, 3f12.5)
99991 format (//2x, 6a4, 14hexceeded limit, i3, 14h on number of ,
     * 9hknots at , e18.8, 13h with error =, e14.4)
      end
