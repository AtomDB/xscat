#include <stdio.h>
#include <math.h>

#include "messages.h"
#include "fitsio.h"
#include "xscat.h"

void read_ZDA_fits(char *ZDA_File, int ZDAnum, struct DUSTMODEL **dm, 
		   int *status) {

  char colname[MAXSTRLEN];
  fitsfile *fptr = NULL;
  int anyn, ii;
  int col;
  long nrows;
  float fval;
  struct DUSTMODEL *ptr;

  ptr = *dm;
  ptr->ZDAtype = ZDAnum;
  ZDAnum++;
  ptr->WDtype = -1;

  if (fits_open_file(&fptr, ZDA_File, READONLY,status)) {
    message("read_ZDA_fits","Error opening file %s",ZDA_File);
    fits_report_error(stderr, *status);
  }
    
  if (fits_movnam_hdu(fptr, ANY_HDU, "ZDADAT", 0, status)) {
    message("read_ZDA_fits","No HDU named ZDADAT as expected.");
    fits_report_error(stderr, *status);
  }
  
  if (fits_get_num_rows(fptr,&nrows,status)) {
    message("read_ZDA_fits","Cannot get number of rows in file?");
    fits_report_error(stderr, *status);
  }
  
  fits_get_colnum(fptr,CASEINSEN,"INAME",&col,status);
  fits_read_col(fptr,TINT,col,ZDAnum,1,1,0,&(ptr->ZDA_iName),&anyn,status);

  fits_get_colnum(fptr,CASEINSEN,"ITYPE",&col,status);
  fits_read_col(fptr,TINT,col,ZDAnum,1,1,0,&(ptr->ZDA_iType),&anyn,status);

  fits_get_colnum(fptr,CASEINSEN,"A",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_AperH = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"AMIN",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_Amin = (double) fval;
  ptr->Amin = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"AMAX",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_Amax = (double) fval;
  ptr->Amax = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"C0",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_C0 = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"B0",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_B[0] = (double) fval;

  ptr->ZDA_A[0] = 0.0;
  ptr->ZDA_M[0] = 0.0;
  for (ii=1;ii<5;ii++) {
    sprintf(colname, "B%d",ii);
    fits_get_colnum(fptr,CASEINSEN,colname,&col,status);
    fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
    ptr->ZDA_B[ii] = (double) fval; 

    sprintf(colname, "A%d",ii);
    fits_get_colnum(fptr,CASEINSEN,colname,&col,status);
    fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
    ptr->ZDA_A[ii] = (double) fval;

    sprintf(colname, "M%d",ii);
    fits_get_colnum(fptr,CASEINSEN,colname,&col,status);
    fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
    ptr->ZDA_M[ii] = (double) fval;
  }

  fits_get_colnum(fptr,CASEINSEN,"ABPPMC",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_ABppmC = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"ABPPMO",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_ABppmO = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"ABPPMSI",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_ABppmSi = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"ABPPMMG",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_ABppmMg = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"ABPPMFE",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_ABppmFe = (double) fval;

  fits_get_colnum(fptr,CASEINSEN,"ABPPMN",&col,status);
  fits_read_col(fptr,TFLOAT,col,ZDAnum,1,1,0,&fval,&anyn,status);
  ptr->ZDA_ABppmN = (double) fval;

  ptr->abppM = 1.e-6*(ptr->ZDA_ABppmC + ptr->ZDA_ABppmO + ptr->ZDA_ABppmSi + 
		      ptr->ZDA_ABppmMg + ptr->ZDA_ABppmFe + ptr->ZDA_ABppmN);
  ptr->atomweight = 1.e-6*(ptr->ZDA_ABppmC*12.0 + ptr->ZDA_ABppmO*16.0 + 
			   ptr->ZDA_ABppmSi*28.09 + ptr->ZDA_ABppmMg*24.31 + 
			   ptr->ZDA_ABppmFe*55.85 + ptr->ZDA_ABppmN*14.01)/
    ptr->abppM;
  ptr->Z = ptr->abppM * ptr->atomweight/1.4;
  
  if (ptr->ZDA_iType == 1) ptr->rho = 2.24; /* % PAHs; g/cm^3 */
  if (ptr->ZDA_iType == 2) {
    if (ptr->ZDA_iName <=6) ptr->rho = 2.24; /* % Graphite; g/cm^3 */
    if (ptr->ZDA_iName > 6) ptr->rho = 1.83; /* % Amorphous Carbon; g/cm^3 */
  }
  if (ptr->ZDA_iType == 3) ptr->rho = 3.5; /* % Silicates; g/cm^3 */
  if (ptr->ZDA_iType == 4) ptr->rho = 3.5; /* % Silicates; g/cm^3 */
  if (ptr->ZDA_iType == 5) {
    if (ptr->ZDA_iName == 4) ptr->rho = 1.05; /* % COMP-GR-S; g/cm^3 */
    if (ptr->ZDA_iName == 5) ptr->rho = 1.05; /* % COMP-GR-FG; g/cm^3 */
    if (ptr->ZDA_iName == 6) ptr->rho = 1.89; /* % COMP-GR-B; g/cm^3 */
    if (ptr->ZDA_iName == 10)ptr->rho = 0.84; /* % COMP-AC-S; g/cm^3 */
    if (ptr->ZDA_iName == 11)ptr->rho = 1.05; /* % COMP-AC-FG; g/cm^3 */
    if (ptr->ZDA_iName == 12)ptr->rho = 1.05; /* % COMP-AC-B; g/cm^3 */
    if (ptr->ZDA_iName == 13)ptr->rho = 1.05; /* % COMP-NC-S; g/cm^3 */
    if (ptr->ZDA_iName == 14)ptr->rho = 1.05; /* % COMP-NC-FG; g/cm^3 */
    if (ptr->ZDA_iName == 15)ptr->rho = 1.05; /* % COMP-NC-B; g/cm^3 */
  }
  ptr->norm = 1.0;

  if (*status != 0) return;

  /* Close the file */
  fits_close_file(fptr, status);
}
