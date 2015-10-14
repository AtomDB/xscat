      real function interp_function(x)

      implicit none

      integer kdimen, ndimen
      parameter(kdimen=5000000, ndimen=12)

      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     *     errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     *     xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     *     xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     *     rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      common/interp_vals/xknots, coefs
      integer knots, iknot
      integer level, nbreak, kntdim, npardm, ibreak, interp, left
      integer maxaux, maxknt, maxpar, maxstk, npar, nstack, leftx
      integer nintrp
c
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     *     dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
c     kntdim - kdimen, name changed to put in common
c     npardm - ndimen, name changed to put in common
      common /resulz/ error, knots
c                      knots = final no. of knots, includes b as one.
c                      error = estimate of error actually achieved.
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
      common /saveit/ iknot
      real ppoly, x
      external ppoly
      
      interp_function  = ppoly(log10(x), xknots, coefs, kdimen, ndimen)

      return
      end

c**************************************************************************

      subroutine build_interp_function(func, lower, upper, epsilon)

      implicit none

      integer kdimen, ndimen
      parameter(kdimen=5000000, ndimen=12)

      real lower, upper, epsilon
      real a, accur, b, bleft, bright, charf, ddtemp, dsctol, error,
     *     errori, factor, fintrp, fleft, fright, norm, xbreak, xdd,
     *     xintrp, xleft, xright
      dimension xbreak(20), dbreak(20), bleft(20), bright(20)
      dimension xleft(50), xright(50), factor(12), fmesge(6)
      dimension ddtemp(12,12), fintrp(10), fleft(6), fright(6),
     *     xdd(12), xintrp(10)
      integer both, break, dbreak, degree, edist, fmesge, right,
     *     rightx, smooth
      logical discrd, fatal, finish
      real xknots(kdimen), coefs(kdimen,ndimen)
      common/interp_vals/xknots, coefs
      real func
      external func
      integer knots, iknot
      integer level, nbreak, kntdim, npardm, ibreak, interp, left
      integer maxaux, maxknt, maxpar, maxstk, npar, nstack, leftx
      integer nintrp
c
      common /inputz/ a, b, accur, norm, charf, xbreak, bleft, bright,
     *     dbreak, degree, smooth, level, edist, nbreak, kntdim, npardm
c     kntdim - kdimen, name changed to put in common
c     npardm - ndimen, name changed to put in common
      common /resulz/ error, knots
c                      knots = final no. of knots, includes b as one.
c                      error = estimate of error actually achieved.
      common /kontrl/ dsctol, errori, xleft, xright, break, both,
     * factor, fmesge, ibreak, interp, left, maxaux, maxknt, maxpar,
     * maxstk, npar, nstack, right, discrd, fatal, finish
      common /comdif/ ddtemp, fintrp, fleft, fright, xdd, xintrp,
     * leftx, nintrp, rightx
      common /saveit/ iknot

      accur = epsilon
      smooth = 0
      degree = 2
      norm = 2
      level = -1
c      level = 1
      a = log10(lower)
      b = log10(upper)
      charf = (b - a)/500
      call adapt(func, xknots, coefs, kdimen, ndimen)

      return
      end
