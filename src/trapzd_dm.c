#include <math.h>

#include "xscat.h"
#include "messages.h"

double dust_mass_dist(struct DUSTMODEL *dm, double Agr);

double trapzd_dm(struct DUSTMODEL *dm, double a, double b, int n) {
  double x,tnm,sum,del;
  static double s;
  int it,j;

  if (n == 1) {
    return (s=0.5*(b-a)*(dust_mass_dist(dm, a)+dust_mass_dist(dm, b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += dust_mass_dist(dm, x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
