#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "xscat.h"
#include "messages.h"

double dust_size_dist(struct DUSTMODEL *dm, double Agr);
double dust_mass_dist(struct DUSTMODEL *dm, double Agr); 
double ZDA_graindist(struct DUSTMODEL *dm, double a);
double WD_graindist(struct DUSTMODEL *dm, double a);

double dust_size_dist(struct DUSTMODEL *dm, double Agr) {

/*
  Goal is to return the unnormalized dust size distribution, 
  n(a)da with units of grains/H atom with sizes between a, a+da.
  Input is dust grain size, in um; output is grains/um.
*/

  double result = -1.0;

  if (dm->WDtype != -1 ) result = 1.e-4*WD_graindist(dm, Agr/1.e4); 
  if (dm->ZDAtype != -1) result = ZDA_graindist(dm, Agr);

  if ((dm->WDtype == -1)&&(dm->ZDAtype == -1)) {
    if (dm->DustModel==MRN_Slc) result = pow(Agr,-3.5);
    if (dm->DustModel==MRN_Gra) result = pow(Agr,-3.5);

    if ((dm->DustModel==WSD_Slc)||(dm->DustModel==WSD_Gra)) {
      if (Agr < 0.5) { 
	result = pow(Agr,-3.5);
      } else { 
	result = pow(Agr,-4);
      } 
    }  
  }
  if ((result != result)||(result-1.0 == result)) {
    message("dust_size_dist","Something wrong");
  }

  if (result < 0) errmess("dust_size_dist","Something wrong");
  return result;

}


double dust_mass_dist(struct DUSTMODEL *dm, double Agr) {
  
  /* Agr is in um.*/
  double result;

  double dnda = dust_size_dist(dm, Agr);  /* grains/H atom/um */
  double dustmass = (4.*M_PI*dm->rho/3.)*pow(Agr/1e4,3); /* g/grain */

  result = dm->norm*dnda*dustmass;
  return result; /* g/H atom/um */
}


double ZDA_graindist(struct DUSTMODEL *dm, double a) { /*  a */

  /*
    %  From Zubko, Dwek & Arendt 2004, ApJS, 152, 211
    %  Output is dN/dA, in [um^-1 H^-1].  
    Input:   A = grain radius (um)
  */

  double gg;
  double Norm, c0, b0, b1, a1, m1, b2, a2, m2, b3, a3, m3, b4, a4, m4;
  double result;

  /* Get values for this particular index */

  Norm = dm->ZDA_AperH;  

  /* c0= dm->ZDA_C0;   
   b0= dm->ZDA_B[0];
   b1= dm->ZDA_B[1];    a1= dm->ZDA_A[1];   m1= dm->ZDA_M[1];
   b2= dm->ZDA_B[2];    a2= dm->ZDA_A[2];   m2= dm->ZDA_M[2];
   b3= dm->ZDA_B[3];    a3= dm->ZDA_A[3];   m3= dm->ZDA_M[3];
   b4= dm->ZDA_B[4];    a4= dm->ZDA_A[4];   m4= dm->ZDA_M[4]; */


   /* Fixed ZDA.fits to have uniform values. */
  /* gg = c0 + b0*log10(a) - 
     b1 * pow( fabs(log10(a/a1)), m1) - 
     b2 * pow( fabs(log10(a/a2)), m2) -
     b3 * pow( fabs(a-a3), m3) - 
     b4 * pow( fabs(a-a4), m4); */

   gg = dm->ZDA_C0 + dm->ZDA_B[0]*log10(a) - 
     dm->ZDA_B[1] * pow( fabs(log10(a/dm->ZDA_A[1])), dm->ZDA_M[1]) - 
     dm->ZDA_B[2] * pow( fabs(log10(a/dm->ZDA_A[2])), dm->ZDA_M[2]) - 
     dm->ZDA_B[3] * pow( fabs(a-dm->ZDA_A[3]), dm->ZDA_M[3]) - 
     dm->ZDA_B[4] * pow( fabs(a-dm->ZDA_A[4]), dm->ZDA_M[4]);

   /* 
      message("calc_dustmodel","%e %e %e %e %e %e %e %e",
              a1,a2,b1,b2,c0,b0,m1,m2);
   */
   result = Norm*pow(10.0,gg);
   if ((result != result)||(result-1.0 == result)) {
     message("calc_dustmodel","Houston, we have a problem, %e",gg);
   }
   return result;
}

double WD_graindist(struct DUSTMODEL *dm, double a) {
  /* Input: a (in units of cm) 
     Output: (1/nH) dn_gr/da, where n_gr(a) is the number density of 
             grains with radius < a and n_H is the H nucleus number density 
             units of cm^-1 */

  double alphagarr[] = 
    {-2.25,-2.17,-2.04,-1.91,-1.84,-1.72,-1.54,-2.26,-2.16,-2.01,
       -1.83,-1.64,-2.35,-2.12,-1.94,-1.61,-2.62,-2.52,-2.36,-2.09,
     -1.96,-2.80,-2.67,-2.45,-1.90,-2.91,-2.99,4.43,-2.94,-2.82,4.16,-2.79};
    double betagarr[] = 
      {-0.0648,-0.0382,-0.111,-0.125,-0.132,-0.322,-0.165,-0.199,-0.0862,
       -0.0973,-0.175,-0.247,-0.668,-0.67,-0.853,-0.722,-0.0144,-0.0541,
       -0.0957,-0.193,-0.813,0.0356,0.0129,-0.00132,-0.0517,
       0.895,2.46,0.0,5.22,9.01,0.0,1.12};
    double atgarr[] = 
      {0.00745,0.00373,0.00828,0.00837,0.00898,0.0254,.0107,0.0241,
       0.00867,0.00811,0.0117,0.0152,0.148,0.0686,0.0786,0.0418,0.0187,
       0.0366,0.0305,0.0199,0.0693,0.0203,0.0134,0.0275,0.012,
       0.578,0.0980,0.00322,0.373,0.392,0.342,0.0190};
    double acgarr[] = 
      {0.606,0.586,0.543,0.499,0.489,0.438,0.428,0.861,0.803,0.696,
       0.604,0.536,1.96,1.35,0.921,0.72,5.74,6.65,6.44,4.6,3.48,3.43,
       3.44,5.14,7.28,
       1.21,0.641,0.285,0.349,0.269,0.0493,0.522};
    double cgarr[] = 
      {9.94e-11,3.79e-10,5.57e-11,4.15e-11,2.90e-11,3.20e-12,9.99e-12,
       5.47e-12,4.58e-11,3.96e-11,1.42e-11,5.83e-12,4.82e-14,3.65e-13,
       2.57e-13,7.58e-13,6.46e-12,1.08e-12,1.62e-12,4.21e-12,2.95e-13,
       2.74e-12,7.25e-12,8.79e-13,2.86e-12,
       7.12e-17,3.51e-15,9.57e-24,9.92e-17,6.20e-17,3.05e-15,8.36e-14};
    double alphasarr[] = 
      {-1.48,-1.46,-1.43,-1.41,-2.1,-2.1,-2.21,-2.03,-2.05,-2.06,-2.08,
       -2.09,-1.57,-1.57,-1.55,-1.59,-2.01,-2.11,-2.05,-2.1,-2.11,-1.09,
       -1.14,-1.08,-1.13,-2.45,-2.49,-2.70,-2.34,-2.36,-2.44,-2.26};
    double betasarr[] = 
      {-9.34,-10.3,-11.7,-11.5,-0.114,-0.0407,0.3,0.668,0.832,0.995,
       1.29,1.58,1.1,1.25,1.33,2.12,0.894,1.58,1.19,1.64,2.1,-0.37,
       -0.195,-0.336,-0.109,0.125,0.345,2.18,-0.243,-0.113,0.254,-3.46};
    double atsarr[] = 
      {0.172,0.174,0.173,0.171,0.169,0.166,0.164,0.189,0.188,0.185,
       0.184,0.183,0.198,0.197,0.195,0.193,0.198,0.197,0.197,0.198,0.198,
       0.218,0.216,0.216,0.211,0.191,0.184,0.198,0.184,0.182,0.188,0.216};
    double csarr[] = 
      {1.02e-12,1.09e-12,1.27e-12,1.33e-12,1.26e-13,1.27e-13,1.e-13,
       5.2e-14,4.81e-14,4.7e-14,4.26e-14,3.94e-14,4.24e-14,4.e-14,
       4.05e-14,3.2e-14,4.95e-14,3.69e-14,4.37e-14,3.63e-14,3.13e-14,
       1.17e-13,1.05e-13,1.17e-13,1.04e-13,
       1.84e-14,1.78e-14,7.29e-15,3.18e-14,3.03e-14,2.24e-14,3.16e-14};
    double bc5arr[] = 
      {0.,1.,2.,3.,4.,5.,6.,0.,1.,2.,3.,4.,0.,1.,2.,3.,0.,1.,2.,3.,4.,
       0.,1.,2.,3.,0.0,1.0,2.0,0.0,0.5,1.0,0.0};

    double alphag = alphagarr[dm->WDtype/2];
    double betag  = betagarr[dm->WDtype/2];
    double atg    = atgarr[dm->WDtype/2]*1.e-4; /* convert to cm */
    double acg    = acgarr[dm->WDtype/2]*1.e-4; /* convert to cm */
    double cg     = cgarr[dm->WDtype/2];
    double alphas = alphasarr[dm->WDtype/2];
    double betas  = betasarr[dm->WDtype/2];
    double ats    = atsarr[dm->WDtype/2]*1.e-4; /* convert to cm */
    double acs    = 1.e-5;
    double cs     = csarr[dm->WDtype/2];
    double bc5    = bc5arr[dm->WDtype/2];
    double dnda;

    /* Calculate sizes */

    if (dm->WDtype/2.0 == dm->WDtype/2) { /* Silicate */
      dnda=(cs/a)*pow(a/ats,alphas);
      if (betas >= 0.) {
	dnda=dnda*(1.+betas*a/ats);
      } else { 
	dnda=dnda/(1.-betas*a/ats);
      }
      if (a > ats) dnda=dnda*exp( pow((ats-a)/acs,3.));
    } else { /* Graphite */
      dnda=(cg/a)*pow(a/atg,alphag);
      if (betag >= 0.) {
	dnda=dnda*(1.+betag*a/atg);
      } else {
	dnda=dnda/(1.-betag*a/atg);
      }
      if (a > atg) dnda=dnda*exp( pow((atg-a)/acg,3.));

      /* Add in PAHs if they are more than 1e-4 of the total. */
      double a01=3.5e-8;
      double a02=3.e-7;
      double sig=0.4;
      double b1=2.0496e-7;
      double b2=9.6005e-11;
      double dndavsg= (b1/a) * exp(-0.5* pow(log(a/a01)/sig,2) )+
	(b2/a)*exp(-0.5* pow(log(a/a02)/sig,2));
      if (dndavsg >= 0.0001*dnda) dnda=dnda+bc5*dndavsg;
    }

    return dnda;
}
