#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "xscat.h"
#include "messages.h"

void expand_env(char *string);
void read_ZDA_fits(char *ZDA_File, int ZDAnum, struct DUSTMODEL **dm, 
		   int *status);
double qtrap_dm(struct DUSTMODEL *dm, double a, double b);

void load_dustmodel(struct PARAMETERS *params, int *status) {

  int indx, dustnum, iDust;
  struct DUSTMODEL *dm;
  double TotMass;
  char tmpstr[MAXSTRLEN];

  for (iDust=0; iDust<params->DustModelNComp; iDust++) {
    dustnum = params->DustModelArr[iDust];

    params->dm[iDust] = (struct DUSTMODEL *) malloc(sizeof(struct DUSTMODEL));
    dm = params->dm[iDust];

    dm->DustModel = dustnum;
    dm->norm = 1.0;
    if ((dustnum == MRN_Slc)||(dustnum == WSD_Slc)) {
      dm->ZDAtype = -1;
      dm->WDtype = -1;
      strcpy(dm->name,"Silicate");
      dm->OptConType = 3;
      dm->abppM = 33.e-6;
      dm->Amin = 0.0050; /* um */
      if (dustnum == MRN_Slc) dm->Amax = 0.250;  /* um */
      if (dustnum == WSD_Slc) dm->Amax = 2.000;  /* um */
      dm->atomweight = 172.0;
      dm->Z = dm->abppM * dm->atomweight/1.4;
      dm->rho = 3.3; /* g/cm^3 */
      /*
	Now to calculate normalization factor.  Should be equal to 
	Z*1.4/(TotMass*N_A), where TotMass is the integral over the dust mass 
	distribution, N_A is Avagadro's number.  
      */
      TotMass = qtrap_dm(dm, dm->Amin,dm->Amax); 
      dm->norm = dm->Z*1.4/(TotMass*6.023e23); /* unitless */
    }

    if ((dustnum == MRN_Gra)||(dustnum == WSD_Gra)) {
      strcpy(dm->name,"Graphite");
      dm->OptConType = 2;
      dm->abppM = 270.e-6;
      dm->Amin = 0.0050; /* um */
      if (dustnum == MRN_Gra) dm->Amax = 0.250;  /* um */
      if (dustnum == WSD_Gra) dm->Amax = 2.000;  /* um */
      dm->atomweight = 12.0;
      dm->Z = dm->abppM * dm->atomweight/1.4;
      dm->rho = 2.2; /* g/cm^3 */
      dm->ZDAtype = -1;
      dm->WDtype = -1;
      
      TotMass = qtrap_dm(dm, dm->Amin,dm->Amax); 
      dm->norm = dm->Z*1.4/(TotMass*6.023e23); /* unitless */
    }

    /* Zubko, Dwek & Arendt Models*/
    if ((dustnum >= ZDAmin)&&(dustnum <= ZDAmax)) {
      strncpy(tmpstr,"$XSCAT/inputs/ZDA.fits",MAXSTRLEN);
      expand_env(tmpstr);
      dm->norm = 1.0;
      read_ZDA_fits(tmpstr, dustnum-ZDAmin, &dm, status);
      if (dm->ZDA_iType == 1) dm->OptConType = 6; /* PAH */
      if (dm->ZDA_iType == 2) dm->OptConType = 4; /* Graphite/Amorph */
      if (dm->ZDA_iType == 3) dm->OptConType = 3; /* Silicate */
      if (dm->ZDA_iType == 4) dm->OptConType = 3; /* Also Silicate */
      if (dm->ZDA_iType == 5) {                   /* Composite */
      	if (dm->ZDA_iName == 1) dm->OptConType =  0; /* BGS; n/a */
      	if (dm->ZDA_iName == 2) dm->OptConType =  0; /* BGF; n/a */
      	if (dm->ZDA_iName == 3) dm->OptConType =  0; /* BGB; n/a */
      	if (dm->ZDA_iName == 4) dm->OptConType =  7; /* CGS */
      	if (dm->ZDA_iName == 5) dm->OptConType =  8; /* CGF */
      	if (dm->ZDA_iName == 6) dm->OptConType =  9; /* CGB */
      	if (dm->ZDA_iName == 7) dm->OptConType =  0; /* BAS; n/a */
      	if (dm->ZDA_iName == 8) dm->OptConType =  0; /* BAF; n/a */
      	if (dm->ZDA_iName == 9) dm->OptConType =  0; /* BAB; n/a */
      	if (dm->ZDA_iName ==10) dm->OptConType = 10; /* CAS */
      	if (dm->ZDA_iName ==11) dm->OptConType = 11; /* CAF */
      	if (dm->ZDA_iName ==12) dm->OptConType = 12; /* CAB */
      	if (dm->ZDA_iName ==13) dm->OptConType = 13; /* CNS */
      	if (dm->ZDA_iName ==14) dm->OptConType = 14; /* CNF */
      	if (dm->ZDA_iName ==15) dm->OptConType = 15; /* CNB */
      }
    }

    /* Weingartner & Draine Models*/
    if ((dustnum >= WDmin)&&(dustnum <= WDmax)) {
      indx = dustnum - WDmin;
      dm->ZDAtype = -1;
      dm->WDtype = indx; 
      dm->norm = 1.0;
      
      if (indx/2.0 == indx/2) {
	strcpy(dm->name,"Silicate");
	dm->OptConType = 3;  /* Use the ZDA04 Silicates; smoother extension */
	dm->abppM = 36.3e-6;
	/* For LMC models, reduce abundance by 1.6 (WD01 page 304) */
	if (indx > 48 && indx <= 61 ) dm->abppM /= 1.6;
	/* For SMC models, reduce abundance by 4 (WD01 page 304) */
	if (indx > 61 && indx < 64 ) dm->abppM /= 4.0;
	dm->Amin = 0.0003; /* um */
	dm->Amax = 0.5;  /* um */
	dm->atomweight = 172.0;
	dm->Z = dm->abppM * dm->atomweight/1.4;
	dm->rho = 3.5; /* g/cm^3 */
      } else {
	strcpy(dm->name,"Graphite");
	dm->OptConType = 2;
	dm->abppM = 330e-6 * 0.7;
	/* For LMC models, reduce abundance by 1.6 (WD01 page 304) */
	if (indx > 48 && indx <= 61 ) dm->abppM /= 1.6;
	/* For SMC models, reduce abundance by 4 (WD01 page 304) */
	if (indx > 61 && indx < 64 ) dm->abppM /= 4.0;
	dm->Amin = 0.0003; /* um */
	dm->Amax = 1.250;  /* um */
	if (indx==51) dm->Amax = 2.5;  /* um */
	dm->atomweight = 12.0;
	dm->Z = dm->abppM * dm->atomweight/1.4;
	dm->rho = 2.24; /* g/cm^3 */
      }
      TotMass = qtrap_dm(dm, dm->Amin,dm->Amax); 
      dm->norm = dm->Z*1.4/(TotMass*6.023e23); /* unitless */
    }
    
    if (params->Drude) dm->OptConType = 0; /* Use Drude regardless */
  }
}
