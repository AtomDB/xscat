#include <stdlib.h>
#include <stdio.h>

#include "ape/ape_error.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_trad.h"

#include "xscat.h"
#include "messages.h"

enum {
  PF_UNKNOWN_TYPE,
  PF_BOOLEAN_TYPE,
  PF_DOUBLE_TYPE,
  PF_FILE_TYPE,
  PF_FLOAT_TYPE,
  PF_INTEGER_TYPE,
  PF_STRING_TYPE
};

#include <assert.h>

#define TRUE 1
#define FALSE 0

static char *OutputFileName;
static char *DustModelName;
static int DustModel;
static double Emin;
static double dE;
static int NumE;
static double ExtractRadius;
static double Xpos;
static int Interp;
static double Epsilon;
static int Drude;
static int clobber;

static Param_Table_Type Control_Parm_Table [] = {
  {"OutputFileName",   PF_FILE_TYPE,    &OutputFileName},
  {"DustModelName",    PF_STRING_TYPE,  &DustModelName},
  {"DustModel",        PF_INTEGER_TYPE, &DustModel},
  {"Emin",             PF_DOUBLE_TYPE,  &Emin},
  {"deltaE",           PF_DOUBLE_TYPE,  &dE},
  {"NumberOfEnergies", PF_INTEGER_TYPE, &NumE},
  {"ExtractRadius",    PF_DOUBLE_TYPE,  &ExtractRadius},
  {"DustPosition",     PF_DOUBLE_TYPE,  &Xpos},
  {"Interpolate",      PF_BOOLEAN_TYPE, &Interp},
  {"Epsilon",          PF_DOUBLE_TYPE,  &Epsilon},
  {"Drude",            PF_BOOLEAN_TYPE, &Drude},
  {"clobber",          PF_BOOLEAN_TYPE, &clobber},
  {NULL, 0, NULL} };

int strip_strcmp(const char *s1, const char *s2);
void expand_env(char *string);
static int get_parameters(Param_Table_Type *params);

void get_input_data(int argc, char *argv[], Param_File_Type *p_input,
		    struct PARAMETERS *params, int *status) {

  char tmpstr[MAXSTRLEN];
  int ii;

  const char *WDnames[NumWD01] = 
    {"WD3100AGAL","WD3110AGAL","WD3120AGAL","WD3130AGAL","WD3140AGAL",
     "WD3150AGAL","WD3160AGAL","WD4000AGAL","WD4010AGAL","WD4020AGAL",
     "WD4030AGAL","WD4040AGAL","WD5500AGAL","WD5510AGAL","WD5520AGAL",
     "WD5530AGAL","WD4000BGAL","WD4010BGAL","WD4020BGAL","WD4030BGAL",
     "WD4040BGAL","WD5500BGAL","WD5510BGAL","WD5520BGAL","WD5530BGAL",
     "WD2600ALMC","WD2610ALMC","WD2620ALMC","WD2600BLMC","WD2605BLMC",
     "WD2610BLMC","WD2900ASMC"};

  const char *ZDAnames[NumZDA] = 
    {"ZDABGS","ZDABGF","ZDABGB","ZDACGS","ZDACGF",
     "ZDACGB","ZDABAS","ZDABAF","ZDABAB","ZDACAS",
     "ZDACAF","ZDACAB","ZDACNS","ZDACNF","ZDACNB"};

  /* Read input parameters and populate params and debug structures. */
  xscat_par_init(argc, argv, p_input, params);

  strncpy(tmpstr,OutputFileName,MAXSTRLEN);
  expand_env(tmpstr);
  sprintf(params->OutputFileName,"%s.dat",tmpstr);

  strncpy(params->DustModelName,DustModelName,MAXSTRLEN);

  params->DustModel  = DustModel;
  params->Emin       = Emin;
  params->dE         = dE;
  params->NumE       = NumE;
  params->ExtractRadius = ExtractRadius;
  params->Xpos       = Xpos;
  params->Interp     = Interp;
  params->Epsilon    = Epsilon;
  params->Drude      = Drude;
  params->clobber    = clobber;

  /* Now parse out which dust model we're playing with here */
  if (DustModel >= 0) {
    params->DustModelNComp = 1;
    params->DustModelArr[0] = DustModel;
  } else {
    if (strcmp(DustModelName,"MRN") == 0) {
      params->DustModelNComp = 2;
      params->DustModelArr[0] = 0;
      params->DustModelArr[1] = 1;
    }
    if (strcmp(DustModelName,"WD01") == 0) { /* baseline WD=WD3160AGAL */
      params->DustModelNComp = 2;
      params->DustModelArr[0] = 6 + 0 + WDmin;
      params->DustModelArr[1] = 6 + 1 + WDmin;
    }
    if (strcmp(DustModelName,"WSD") == 0) {
      params->DustModelNComp = 2;
      params->DustModelArr[0] = 4;
      params->DustModelArr[1] = 5;
    }
    
    for (ii=0; ii<NumZDA; ii++) {
      if (strcmp(DustModelName,ZDAnames[ii]) == 0) {
	params->DustModelNComp = 5;
	params->DustModelArr[0] = 5*ii + 0 + ZDAmin;
	params->DustModelArr[1] = 5*ii + 1 + ZDAmin;
	params->DustModelArr[2] = 5*ii + 2 + ZDAmin;
	params->DustModelArr[3] = 5*ii + 3 + ZDAmin;
	params->DustModelArr[4] = 5*ii + 4 + ZDAmin;
      }
    }
      
    for (ii=0;ii<NumWD01;ii++) {
      if (strcmp(DustModelName,WDnames[ii]) == 0) {
	params->DustModelNComp = 2;
	params->DustModelArr[0] = 2*ii + 0 + WDmin;
	params->DustModelArr[1] = 2*ii + 1 + WDmin;
      }
    }
  }
}

int strip_strcmp(const char *s1, const char *s2) {
  
  /* Strips any whitespace from string s1, and compares it to s2. */

  char strip_string[MAXSTRLEN];

  /*  strcpy(strip_string,s1); */
  sscanf(s1,"%s",strip_string);
  return strcmp(strip_string,s2);
}

static int get_parameters(Param_Table_Type *params) {

  int status=eOK;

  for (; NULL!=params->name && eOK == status; ++params) {
    /* Handle strings a little carefully; ape allocates just enough space. */

    char * tmp_str = 0;
    if (0 != tmp_str) fprintf(stderr, "You're wrong!!!!\n");

    switch (params->type) {
      case PF_BOOLEAN_TYPE: {
        int *value = (int *) params->value;
        char bool_value='\0';
        status = ape_trad_query_bool(params->name, &bool_value);
        *value = bool_value;
        break;
      }
      case PF_DOUBLE_TYPE:
        status = ape_trad_query_double(params->name, (double *) params->value);
        break;
      case PF_FILE_TYPE:
        status = ape_trad_query_file_name(params->name, &tmp_str);
        break;
      case PF_FLOAT_TYPE:
        status = ape_trad_query_float(params->name, (float *) params->value);
        break;
      case PF_INTEGER_TYPE:
        status = ape_trad_query_int(params->name, (int *) params->value);
        break;
      case PF_STRING_TYPE:
        /* Try to get string, matching the case of enumerated parameters. */
        status = ape_trad_query_string_case(params->name, &tmp_str, eEnumCase);
        /* If the parameter wasn't enumerated, that's OK, just ignore the error code. */
        if (eRangeNoEnum == status) status = eOK;
        break;
      default:
        status = -1;
    }
    /* If tmp_str was used, copy it over to the actual parameter value. */
    if (eOK == status && 0 != tmp_str) {
      char * real_str = calloc(MAXSTRLEN, sizeof(char));
      if (0 == real_str) status = eDynAllocFailed;
      if (eOK == status) {
        strncpy(real_str, tmp_str, MAXSTRLEN - 1);
        *(char **) params->value = real_str;
      }
      free(tmp_str);
    }
  }
  return status;
}



void xscat_par_init(int argc, char *argv[], Param_File_Type *p_input,
		     struct PARAMETERS *params) {

  int status=0;

  /* ************************** */
  /* Check the input parameters */
  /* ************************** */

  if (eOK != (status = ape_trad_init (argc, argv))) {
    errmess ("xscat_par_init", "Error opening parameter file.\n");
  }

  if (eOK == status) {
    if (eOK != (status = get_parameters (Control_Parm_Table))) {
      ape_trad_close (0);
      errmess("xscat_par_init", 
	      "Error getting parameters (Ape error code %d).\n", status);
    } else {
      ape_trad_close (0);
    }
  }

}

