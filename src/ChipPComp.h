/* This file is part of ChipPComp,
   http://github.com/alexjgriffith/ChipPComp/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define R_NO_REMAP
#define MAX_BED_BUFFER_SIZE 1024

#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include <Rinternals.h>


void pileup(char ** filename,char ** chro,int *start,
	    int *end,int *length,int *scores);


void unityOutput(int * intChr, int * intSummit, int * name, int * l1,
		 int * l2, int * peaklength, int* peakwidth,
		 int* retChr, int * retSummit, int * retMatrix);

static R_NativePrimitiveArgType pileup_t[]={STRSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType unityOutput_t[]={INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};


static R_CMethodDef cMethods[]={
  {"pileup",(DL_FUNC) &pileup,6, pileup_t},
  {"unityOutput",(DL_FUNC) &unityOutput,10,unityOutput_t},
  {NULL,NULL,0,NULL}
};

void R_init_ChipPComp(DllInfo *info);
