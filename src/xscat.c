#include <stdlib.h>
#include <stdio.h>
/* #include "fitsio.h" */

#include "xscat.h"

void get_input_data(int argc, char *argv[], Param_File_Type *p_input,
		    struct PARAMETERS *params, int *status);
void load_dustmodel(struct PARAMETERS *params, int *status);
void set_dust_parm(struct PARAMETERS *params, int *status);
double qtrap_da(struct PARAMETERS *params, struct DUSTMODEL *dm, double Energy);
double qromb_da(struct PARAMETERS *params, struct DUSTMODEL *dm, double Energy);

int main(int argc, char *argv[]) {

  struct PARAMETERS params;
  struct DUSTMODEL *dm;
  int status=0, iE, iDust;
  FILE *fp;
  double result[MAXDUSTCOMP], totscat;
  double Energy;

  /* Get needed inputs: dust type, output file name, etc */
  get_input_data(argc, argv, NULL, &params, &status);
  if (status != 0) {
    printf("xscat: Error reading input data.  Exiting without running.\n");
    return(1); /* error reading input */
  }

  if (!params.clobber) {
    fp = fopen(params.OutputFileName,"r");
    if (fp != NULL) {
      printf("xscat: File %s exists and clobber=yes. Exiting.\n", 
	     params.OutputFileName);
      return(1); /* error reading input */
    }
  }
  fp = fopen(params.OutputFileName,"w");

    /* Read in files, set things up, etc */ 
  load_dustmodel(&params, &status);
  if (status != 0) {
    printf("xscat: Error doing initialization.  Exiting without running.\n");
    return(1); /* error reading input */
  }

  for (iE=0;iE<params.NumE;iE++) {
    Energy = params.Emin + params.dE*iE;
    fprintf(fp,"%e ",Energy);

    /* Calculate the total cross section for scattering */
    /* Integrate \int \int \int da dx dt n(a) dsigma(a, phi)/(1-x)^2 */
    /* where phi = theta/(1-x) */
    /* Phi max = 5 * 624 arcsec / (E(keV) * (a/0.1 um)) */

    totscat = 0.0;
    for (iDust=0;iDust<params.DustModelNComp;iDust++) {
      dm = params.dm[iDust];
      /* result[iDust] = qromb_da(&params, dm, Energy); */
      result[iDust] = qromb_da(&params, dm, Energy);
      fprintf(fp,"%e ", result[iDust]);
      totscat += result[iDust];
    }
    fprintf(fp,"%e \n", totscat);
  }
  fclose(fp);

  return(status); 
}
