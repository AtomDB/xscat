#include <stdlib.h>
#include <stdio.h>

#include "xscat.h"

void wrap_receive_sigmaint_(float *Energy_fl, float *size_a_fl, 
			    float *size_m_fl, float *xpos_fl, 
			    int *Interp, float *Norm_fl, int *Drude, 
			    int *DustType, int *MantleType, float *rho_fl, 
			    float *Tmax_fl, float *epsilon_fl, 
			    float *result_fl);

void wrap_send_sigmaint(struct PARAMETERS *params, struct DUSTMODEL *dm,
			double Energy, double size_a, double size_m, 
			double thetamax, double *IntenInt) {

  int Interp, Drude, DustType, MantleType;
  float Energy_fl, size_a_fl, size_m_fl, Norm_fl, rho_fl, Tmax_fl;
  float epsilon_fl, xpos_fl, result_fl;

  Energy_fl = (float) Energy;
  size_a_fl = (float) size_a;
  size_m_fl = (float) 0.0;  /* Skip mantles for the moment */
  xpos_fl   = (float) params->Xpos;
  Tmax_fl   = (float) thetamax;
  Interp    = params->Interp;
  Norm_fl   = 1.0;
  Drude     = params->Drude;
  DustType  = dm->OptConType;
  MantleType= 0; 
  rho_fl    = dm->rho; /* Only used, I think, in Drude approx */
  epsilon_fl= (float) params->Epsilon;
  
  wrap_receive_sigmaint_(&Energy_fl, &size_a_fl, &size_m_fl, &xpos_fl, 
			 &Interp, &Norm_fl, &Drude, 
			 &DustType, &MantleType, &rho_fl, 
			 &Tmax_fl, &epsilon_fl, &result_fl);
  *IntenInt = (double) result_fl;

}
