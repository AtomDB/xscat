#include <XSFunctions/Utilities/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"
#include <cmath>
#include <iomanip>

//#include <stlToCArrays.h>
#include <sstream>
#include <iostream>
#include <string>
#include <fitsio.h>
#include <fstream>

using namespace XSutility;
using namespace std;

#include "xscat.h"
void xscat_init(int *Nmodel, struct SIGMA_TABLE *xscatdata[]);
void calc_xscat(const RealArray& energyArray, 
		RealArray& fluxmult, 
		struct SIGMA_TABLE *xscat_dat, 
		Real NH, Real xpos, Real rext);

extern "C" 
void xscatmodel(const RealArray& energyArray, 
	   const RealArray& params, 
	   int spectrumNumber, 
	   RealArray& fluxmult, 
	   RealArray& fluxmultError, 
	   const string& initString) {
  
  /* energy : size Nflux+1 (has low and high bins, in keV)
     parameter: all the parameters
     spectrum: number of component being calculated.  Ignored.
     flux: output flux
     fluxError: ignore
     init: initialization string - may hold location of input file.
  */

  std::ostringstream errMsg;
  static bool isFirst = true;
  static int Nmodel;
  static struct SIGMA_TABLE *xscat_dat[MAXMODELS];
  
  const Real NH   = params[0];
  const Real xpos = params[1];
  const Real rext = params[2];
  int   dmodel    = (int) params[3];  

  if (isFirst) {
    Nmodel = 0;
    xscat_init(&Nmodel, xscat_dat);
    isFirst = false;
  }

  if ((dmodel < 1)||(dmodel > Nmodel)) {
    errMsg << "\nxscat:\nDust model #" << 
      (int) params[5] << " must be between 1 and " << Nmodel << "\n";
    throw YellowAlert(errMsg.str());
  }
  
  size_t N(energyArray.size());
  fluxmult.resize(N-1);

  /* ********************************* */
  /* Now do the scattering calculation */
  /* ********************************* */

  calc_xscat(energyArray, fluxmult, xscat_dat[dmodel-1], NH, xpos, rext);

  /* And we're done */

  fluxmultError.resize(0);
  return;
}


