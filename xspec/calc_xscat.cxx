#include <XSUtil/FunctionUtils/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/FunctionUtils/FunctionUtility.h>
#include <XSUtil/FunctionUtils/xsFortran.h>
#include <XSFunctions/functionMap.h>
/* #include "xsTypes.h" */
#include <cmath>
#include <iomanip>

//#include <stlToCArrays.h>
#include <sstream>
#include <iostream>
#include <string>
#include <fitsio.h>
#include <fstream>

#include "xscat.h"
Real interpol_huntf(int n, Real *x, Real *y, Real z);

void calc_xscat(const RealArray& energyArray, 
		RealArray& fluxmult, 
		struct SIGMA_TABLE *xscat_dat, 
		Real NH, Real xpos, Real rext) {

  std::ostringstream errMsg;
  size_t N(energyArray.size());
  Real Eval[MAXENERGY], Sval[MAXENERGY], en, sigma=0.0;
  size_t ii, jj, kk, iE, Nrext, Nxpos;
  int rl, rh, xl, xh;
  Real xl_wt=0.0, xh_wt=0.0, rl_wt=0.0, rh_wt=0.0;

  rl = -1;
  rh = -1;
  xl = -1;
  xh = -1;

  // printf("Inputs: %e, %e, %e\n", NH, xpos, rext);
  
  Nrext = xscat_dat->Nrext;
  Nxpos = xscat_dat->Nxpos;
  
  // Bound the parameters in our list.  Note this assumes a specific ordering
  // to the FITS file.

  for (jj=0;jj<Nxpos-1;jj++) {
    if ( (xscat_dat->xpos[jj*Nrext] <= xpos) &&
	 (xscat_dat->xpos[(jj+1)*Nrext] > xpos) ) {
      xl = jj;
      xh = jj+1;
      xl_wt = (xscat_dat->xpos[(jj+1)*Nrext] - xpos)/
	(xscat_dat->xpos[(jj+1)*Nrext] - xscat_dat->xpos[jj*Nrext]);
      xh_wt = 1 - xl_wt;
    }
  }
  if ((rext >= 0.0)&&(rext < xscat_dat->rext[0])) rext = xscat_dat->rext[0];
  for (kk=0; kk<xscat_dat->Nrext-1; kk++) {
    if ((xscat_dat->rext[kk] <= rext) &&
	(xscat_dat->rext[kk+1] > rext) ) {
      rl = kk;
      rh = kk+1;
      rl_wt = (xscat_dat->rext[kk+1] - rext)/
	(xscat_dat->rext[kk+1] - xscat_dat->rext[kk]);
      rh_wt = 1-rl_wt;
    }
  }
  // printf("%e, %e: Pos: %d %d %d %d\n", xpos, rext, xl, xh, rl, rh);
  // printf("Wt:  %e %e %e %e\n", xl_wt, xh_wt, rl_wt, rh_wt);

  if ((xl<0)||(xh<0)||(rl<0)||(rh<0)) {
    errMsg << "\nxscat:\nDust position " << xpos << 
      "not in range 0->1 or extraction regions " << rext <<
      " out of bounds.\n";
    throw YellowAlert(errMsg.str());
  }
  
  // Build up the best-fit sigma values as a function of energy
  for (ii=0; ii<xscat_dat->Nenergy; ii++) {
    iE = ii*(xscat_dat->Nxpos*xscat_dat->Nrext);
    Eval[ii] = xscat_dat->energy[iE];
    
    Sval[ii] = 
      xl_wt*rl_wt*xscat_dat->sigma[iE + xl*Nrext + rl] +
      xl_wt*rh_wt*xscat_dat->sigma[iE + xl*Nrext + rh] + 
      xh_wt*rl_wt*xscat_dat->sigma[iE + xh*Nrext + rl] +
      xh_wt*rl_wt*xscat_dat->sigma[iE + xh*Nrext + rh];
    /* printf("Sval: %e %e %e\n", Eval[ii], Sval[ii], 
       xscat_dat->sigma[iE + xl*Nrext + rh]); */
  }
  
  for (ii=0; ii<N-1; ii++) {
    en = 0.5*(energyArray[ii] + energyArray[ii+1]);
    if ((en >= Eval[0])&&(en <= Eval[xscat_dat->Nenergy-1])) {
      sigma = interpol_huntf(xscat_dat->Nenergy, Eval, Sval, en);
    } else {
      if (en < Eval[0]) fluxmult[ii] = Sval[0]; // Simplistic...
      if (en > Eval[xscat_dat->Nenergy-1]) // RG scaling, E^-2.
	sigma = Sval[xscat_dat->Nenergy-1]*
	  pow(en/Eval[xscat_dat->Nenergy-1],-2.0);
    }
    fluxmult[ii] = exp(-1.e22*NH*sigma);
    /* printf("Xscat: %e %e %e\n", en, sigma, fluxmult[ii]); */
  }

}

Real interpol_huntf(int n, Real *x, Real *y, Real z) {

  Real grad,d1,df,f;
  int inc, jl, jm, ju;
  
  inc = (x[n-1] > x[0]) ? TRUE : FALSE;
  if (( inc && (z < x[0] || z > x[n-1])) ||
      (!inc && (z > x[0] || z < x[n-1]))) {
    printf("interpol_huntf: Asking for %e, min is %e, max is %e\n",z,
            x[0],x[n-1]);
    printf("interpol_huntf: Cannot extrapolate\n");
    abort();
  }  
  jl = 0;
  ju = n;
  while (ju - jl > 1) {
    jm = (ju + jl) / 2;
    if ((z > x[jm]) == inc) {
      jl = jm;
    } else {
      ju = jm;
    }
  }
  /* ------     Z is sandwiched between JL and JU ------ */
  if ((x[jl] > 0. && x[ju] > 0.) && (y[jl] > 0. && y[ju] > 0.)) {
    grad = (log10(y[ju]) - log10(y[jl])) /
      (log10(x[ju]) - log10(x[jl]));
    df = grad * (log10(z) - log10(x[jl]));
    d1 = log10(y[jl]) + df;
    f = pow(10., d1);
  } else {
    f = y[jl]+(y[ju]-y[jl])*((z-x[jl])/(x[ju]-x[jl]));
  }
  return f;
}
