#include <math.h>
#define EPS 1.0e-4
#define JMAX 20

#include "xscat.h"
#include "messages.h"

double trapzd_dm(struct DUSTMODEL *dm, double a, double b, int n);

double qtrap_dm(struct DUSTMODEL *dm, double a, double b) {

	int j;
	double s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd_dm(dm,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		if (s == 0.0 && olds == 0.0 && j > 6) return s;
		olds=s;
	}
	errmess("qtrap_dm","Too many steps in routine qtrap");
	return 0.0;
}
#undef EPS
#undef JMAX
