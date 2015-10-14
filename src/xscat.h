/* 
  Header File for program "xscat.c"
  written by Randall K. Smith
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#ifndef xscat_h
#define xscat_h

#ifdef __cplusplus
extern "C" {
#endif
  
#define VERSION "1.0"
  
#define MAXLINELENGTH 132
#ifndef MAXSTRLEN
#define MAXSTRLEN 1024
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define MRN_Slc  0
#define MRN_Gra  1
#define WD_Slc   2  /* R_V=3.1, b_C = 6e-5 */
#define WD_Gra   3  
#define WSD_Slc  4  /* Witt, Smith & Dwek 2001 */
#define WSD_Gra  5  /* Witt, Smith & Dwek 2001 */
#define ZDAmin   6
#define ZDAmax  80
#define SMC_Slc 81 /* Model 32 for SMC; WD2001, ApJ, 548, 296 */
#define SMC_Gra 82 
#define WDmin   83
#define WDmax  147  
#define NumWD01 32
#define NumZDA  15

#define MAXDUSTCOMP 10 

  struct DUSTMODEL {
    int DustModel;
    int OptConType;
    char name[MAXSTRLEN];
    double abppM;
    double Amin;
    double Amax;
    /* double *sizefunc; */
    /* double ap[4]; */
    double Z;
    double rho;
    double norm;
    double atomweight;
    /* int henke_N; */
    /* double *henke_E; */
    /* double *henke_F; */
    int WDtype;
    int ZDAtype;
    int ZDA_iName;
    int ZDA_iType;
    float ZDA_AperH, ZDA_Amin, ZDA_Amax, ZDA_C0;
    float ZDA_B[5], ZDA_A[5],ZDA_M[5];
    float ZDA_ABppmC,ZDA_ABppmO,ZDA_ABppmSi,ZDA_ABppmMg,ZDA_ABppmFe,ZDA_ABppmN;
  };

  /* ****************************************************** */
  /* PARAMETERS contains values that are set at the outset  */
  /* of the run, and do not change afterward.               */
  /* ****************************************************** */
  
  struct PARAMETERS {       /* structure to collect all the plasma parameters */
    char OutputFileName[MAXSTRLEN];
    char DustModelName[MAXSTRLEN];
    double Emin;
    double dE;
    double Energy;
    double ExtractRadius;
    double Xpos;
    double Epsilon;
    int NumE;
    int Interp;
    int DustModel;
    int DustModelNComp;
    int DustModelArr[MAXDUSTCOMP];
    int Drude;
    int ScatteringType;
    struct DUSTMODEL *dm[MAXDUSTCOMP];
    int clobber;
  };
  
  typedef struct Param_Table {
    char * name;
    unsigned int type;
    void * value;
  } Param_Table_Type;
  
  typedef void Param_File_Type;
  
  void xscat_par_init(int argc, char *argv[], Param_File_Type *p_input,
		      struct PARAMETERS *params);
  
#ifdef __cplusplus
}
#endif

#endif
