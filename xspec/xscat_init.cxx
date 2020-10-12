#include <XSFunctions/Utilities/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
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

using namespace std;

void read_xscat(const char *filename, int *Nmodel, 
		struct SIGMA_TABLE *xscat_dat[]);

void xscat_init(int *Nmodel, struct SIGMA_TABLE *xscatdata[]) {

  std::ostringstream errMsg;

  string mrnfile, wd1file, zda1file;
  
  *Nmodel=0;
  
  mrnfile  = DIRECTORY "/xscat_MRN_v1.0.0.fits";
  wd1file  = DIRECTORY "/xscat_WD3100AGAL_v1.0.0.fits";
  zda1file = DIRECTORY "/xscat_ZDABAS_v1.0.0.fits";

  try {
    ifstream file_exists(mrnfile.c_str());
      if (!file_exists) {
	errMsg << "\nxscat:\nError finding MRN scattering file: " 
	       << mrnfile.c_str() << "\n";
	throw YellowAlert(errMsg.str());
      }
      read_xscat(mrnfile.c_str(), Nmodel, xscatdata);
  }
  catch (...) {
    return;
  }

  try {
    ifstream file_exists(wd1file.c_str());
      if (!file_exists) {
	errMsg << "\nxscat:\nError finding WD RV=3.1, bC=0.0 scattering file.\n";
	throw YellowAlert(errMsg.str());
      }
      read_xscat(wd1file.c_str(), Nmodel, xscatdata);
  }
  catch (...) {
    return;
  }

  try {
    ifstream file_exists(zda1file.c_str());
      if (!file_exists) {
	errMsg << "\nxscat:\nError finding ZDABAS scattering file.\n";
	throw YellowAlert(errMsg.str());
      }
      read_xscat(zda1file.c_str(), Nmodel, xscatdata);
  }
  catch (...) {
    return;
  }
}

void read_xscat(const char *filename, int *Nmodel, 
		struct SIGMA_TABLE *xscat_dat[]) {
  
  fitsfile *fptr;
  int hdunum, casesen;
  
  struct SIGMA_TABLE *xscat;
  int colnum_energy, colnum_dustpos, colnum_radextr, colnum_sigma;
  char column[MAXSTRLEN],comment[MAXSTRLEN];

  int status, hdutype, anynul;
  int nrows, nxpos, nrext, nenergy;
  int row;
  int iX, iE, iR;
  char EnergyStr[]  = "Energy";
  char DustPosStr[] = "DustPos";
  char RadExtrStr[] = "RadExtr";
  char SigmaStr[]   = "Sigma";
  
  xscat = (struct SIGMA_TABLE *) malloc(sizeof(struct SIGMA_TABLE));

  std::ostringstream errMsg;
  status = 0;

  if (fits_open_file(&fptr, filename, READONLY, &status)) {
    errMsg << status << "\nFailed to open " << filename << "\n";
    throw YellowAlert(errMsg.str());
  }

  hdunum = 2;
  if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) {
    errMsg << status << "\nread_xscat\nError finding second HDU\n";
    throw YellowAlert(errMsg.str());
  }
    
  casesen = FALSE;

  fits_get_colname(fptr, casesen, EnergyStr, column, &colnum_energy, &status);
  fits_get_colname(fptr, casesen, DustPosStr,column, &colnum_dustpos, &status);
  fits_get_colname(fptr, casesen, RadExtrStr, column,&colnum_radextr, &status);
  fits_get_colname(fptr, casesen, SigmaStr, column,  &colnum_sigma, &status);

  if (status) {
    errMsg << status << "\nread_xscat\nError getting column names\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_key(fptr, TINT, "NAXIS2", &nrows, comment, &status)) {
    errMsg << status << "\nread_xscat\nError getting number of rows\n";
    throw YellowAlert(errMsg.str());
  }

  if (fits_read_key(fptr, TINT, "NENERGY", &nenergy, comment, &status)) {
    errMsg << status << "\nread_xscat\nError getting number of rows\n";
    throw YellowAlert(errMsg.str());
  }

  if (fits_read_key(fptr, TINT, "NXPOS", &nxpos, comment, &status)) {
    errMsg << status << "\nread_xscat\nError getting number of rows\n";
    throw YellowAlert(errMsg.str());
  }

  if (fits_read_key(fptr, TINT, "NREXT", &nrext, comment, &status)) {
    errMsg << status << "\nread_xscat\nError getting number of rows\n";
    throw YellowAlert(errMsg.str());
  }

  xscat->Nrows = (size_t) nrows;
  xscat->Nenergy = (size_t) nenergy;
  xscat->Nxpos = (size_t) nxpos;
  xscat->Nrext = (size_t) nrext;
  xscat->energy = (float *) malloc(nrows*sizeof(float));
  xscat->xpos   = (float *) malloc(nrows*sizeof(float));
  xscat->rext   = (float *) malloc(nrows*sizeof(float));
  xscat->sigma  = (float *) malloc(nrows*sizeof(float));

  if (nenergy*nxpos*nrext != nrows) {
    errMsg << status << "\nread_xscat\nNumber of rows" << nrows << 
      "!= Number of energies * positions * extraction radii.  Odd." <<")\n";
    throw YellowAlert(errMsg.str());
  }

  row = 0;
  for (iE=0;iE<xscat->Nenergy;iE++) {
    for (iX=0;iX<xscat->Nxpos;iX++) {
      for (iR=0;iR<xscat->Nrext;iR++) {
	row++;

	if (fits_read_col(fptr, TFLOAT, colnum_energy, row, 1, 1, NULL,
			  &(xscat->energy[row-1]), &anynul, &status)) {
	  errMsg << status << "\nread_xscat\nError reading energy\n";
	  throw YellowAlert(errMsg.str());
	}
      
	if (fits_read_col(fptr, TFLOAT, colnum_dustpos, row, 1, 1, NULL,
			  &(xscat->xpos[row-1]), &anynul, &status)) {
	  errMsg << status << "\nread_xscat\nError reading position values\n";
	  throw YellowAlert(errMsg.str());
	}

	if (fits_read_col(fptr, TFLOAT, colnum_radextr, row, 1, 1, NULL,
			  &(xscat->rext[row-1]), &anynul, &status)) {
	  errMsg << status << 
	    "\nread_xscat\nError reading extraction radius values\n";
	  throw YellowAlert(errMsg.str());
	}
	
	if (fits_read_col(fptr, TFLOAT, colnum_sigma, row, 1, 1, NULL,
			  &(xscat->sigma[row-1]), &anynul, &status)) {
	  errMsg << status << "\nread_xscat\nError reading sigma values\n";
	  throw YellowAlert(errMsg.str());
	}
      }
    }      
  }

  if (fits_close_file(fptr, &status)) {
    errMsg << status << "\nread_xscat\nError closing fits file\n";
    throw YellowAlert(errMsg.str());
  }

  xscat_dat[*Nmodel] = xscat;
  (*Nmodel)++;
}
