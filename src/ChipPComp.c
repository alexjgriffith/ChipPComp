/* This file is part of ChipPComp,
   http://github.com/alexjgriffith/ChipPComp/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include "ChipPComp.h"

// #define _VERBOSE 1
// #define _test 1

void unityOutput(int * intChr, int * intSummit, int * name, int * l1,
                 int * l2, int * peaklength, int* peakwidth, int* retChr,
                 int * retSummit, int * retMatrix)
{
   int i,j,k;
   int nextSummit,temp;
   for(i=0;i<(*peaklength);i++)
     {       
       retChr[i]=intChr[l1[i]-1];
       nextSummit=0;
       j=l2[i]-l1[i]+1;
       if(j==0)	
	 Rf_error("Divide by zero in file<-region.c function<-unityOutput\n");
      
       for(k=(l1[i]-1);k<l2[i];k++)
	 {
	   retMatrix[i* (*peakwidth)+name[k]-1]=1;
	   temp=nextSummit+intSummit[k]/j;
	   nextSummit=temp;
	 }	  	
       retSummit[i]=nextSummit;      
     }  
}

void pileup(char ** filename,char ** chro,int *start,
	    int *end,int *length,int *scores)
{
#ifdef _VERBOSE
  Rprintf("N=2\tfilename=%s\n",*filename);
#endif
#ifndef _test
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  char string[1024];
  int inStart,inEnd;
  int i=0;
  
  int chrC;
  int ed1C;
  int ed2C;

  while(fgets(buffer,1024, f))
    {            
      sscanf(buffer,"%s\t%d\t%d", string,&inStart,&inEnd);
      chrC=strcmp(chro[i],string); // < 0 next peak
                                   // == 0 next compare
                                   // > 0 next reed
      while(chrC<=0){
	ed1C=inStart-end[i]; //  =< 0 -- inStart < end
	                           //  > 0 next read
	ed2C=start[i]-inEnd; // =< 0-- inEnd > start
                                    // > 0 -- next peak
	if(chrC==0){
	  if(ed1C<=0 && ed2C<=0){
	    scores[i]++;
	    break;
	  }
	  else if (ed1C>0){
	    i++;

	  }
	  else{
	    break;
	  }
	    
	}
	else{

	  i++;
	}
	if(i>=*length)
	  break;	
	chrC=strcmp(chro[i],string);
	}
      if(i>=*length)
	break;	
    }

  fclose(f);
#else
  Rprintf("N=2\tfilename=%s\n",*filename);
  Rprintf("N=2\ttesting\n",*filename);
  Rprintf("N=2\t%d\n",*length);
  Rprintf("N=2\t%s\n",chro[0]);
  int i;
  for(i=0;i<*length;i++)
    Rprintf("%s\n",chro[i]);

#endif
}

void R_init_ChipPComp(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
