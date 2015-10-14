#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "xscat.h"
#include "messages.h"

#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

double dust_size_dist(struct DUSTMODEL *dm, double Agr);
void wrap_send_sigmaint(struct PARAMETERS *params, struct DUSTMODEL *dm,
			double Energy, double size_a, double size_m, 
			double thetamax, double *IntenInt);
void polint_da(double xa[], double ya[], int n,double x,double *y,double *dy);
double trapzd_da(struct PARAMETERS *params, struct DUSTMODEL *dm, 
		 double Energy, double a, double b, int n);

double qtrap_da(struct PARAMETERS *params, struct DUSTMODEL *dm, double Energy){

  int j;
  double s,olds;

  olds = -1.0e30;
  for (j=1;j<=JMAX;j++) {
    s=trapzd_da(params, dm, Energy, dm->Amin, dm->Amax, j);
    if (fabs(s-olds) < params->Epsilon*fabs(olds)) return s;
    if (s == 0.0 && olds == 0.0 && j > 6) return s;
    message("qtrap","%e  %e %e %e",s,olds, fabs(s-olds)/fabs(olds),
	    params->Epsilon); 
    olds=s;
  }
  message("quadrature_dustsize","Too many steps in routine qromb");
  return 0.0;
}

double qromb_da(struct PARAMETERS *params, struct DUSTMODEL *dm, double Energy){

  double a, b;
  double ss,dss;
  double s[JMAXP],h[JMAXP+1];
  int j;

  a = dm->Amin;
  b = dm->Amax;

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd_da(params, dm, Energy, a, b, j);
    if (j >= K) {
      polint_da(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= params->Epsilon*fabs(ss)) return ss;
    }
    if (ss!=ss) message("Qd","Qromb: %e,%d,%e,%e",s[j],j,ss,dss);
    /* message("qromb","%e  %e",ss,fabs(dss));  */
    h[j+1]=0.25*h[j];
  }
  message("quadrature_dustsize","Too many steps in routine qromb");
  return 0.0;
}

void polint_da(double xa[], double ya[], int n,double x,double *y,double *dy){

  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;
  
  dif=fabs(x-xa[1]);
  c = (double *) malloc((n+1)*sizeof(double));
  d = (double *) malloc((n+1)*sizeof(double));

  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) message("polint_da","Error in routine polint");
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(c);
  free(d);
}

double trapzd_da(struct PARAMETERS *params, struct DUSTMODEL *dm, 
		 double Energy, double a, double b, int n) {
  double x,tnm,sum,del;
  double FUNCa, FUNCb, FUNCx;
  double inten_a, inten_b, inten_x;
  double thetamax;
  static double s;
  int it,j;

  thetamax = params->ExtractRadius;

  if (n == 1) {
    FUNCa = dust_size_dist(dm, a);
    if (FUNCa > 0.0) {
      wrap_send_sigmaint(params, dm, Energy, a, 0.0, thetamax, &inten_a);
      FUNCa *= inten_a*dm->norm;
    }
    FUNCb = dust_size_dist(dm, b); 
    if (FUNCb > 0.0) {
      wrap_send_sigmaint(params, dm, Energy, b, 0.0, thetamax, &inten_b);
      FUNCb *= inten_b*dm->norm;
    }
    return (s=0.5*(b-a)*(FUNCa+FUNCb));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) {
      FUNCx = dust_size_dist(dm, x);
      if (FUNCx > 0.0) {
	wrap_send_sigmaint(params, dm, Energy, x, 0.0, thetamax, &inten_x);
	if (inten_x != inten_x) message("Qd","Failed in sigma, e=%e",Energy);
	FUNCx *= inten_x*dm->norm;
      }
      sum += FUNCx;
    }
    s=0.5*(s+(b-a)*sum/tnm);
    if (s!=s) message("Qd","How possible? %e, %e, %e, %d",s,a,b,tnm);
    return s;
  }
}

#undef JMAX
#undef JMAXP
#undef K
