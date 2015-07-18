/*************************************************************************

   Program:    countcg
   File:       countcg.c
   
   Version:    V1.1
   Date:       27.09.94
   Function:   Count the conformations in a CG file (ignore reference set)
   
   Copyright:  (c) Dr. Andrew C. R. Martin, 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling,
               University College London,
   Phone:      HOME: +44 (0)372 275775
   EMail:      JANET:    martin@uk.ac.ucl.bioc.bsm
               INTERNET: martin@bsm.bioc.ucl.ac.uk
                         amartin@scitec.adsp.sub.org
               UUCP:     ...cbmehq!cbmuk!scitec!amartin
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and 
   that the source code is provided with the program.

**************************************************************************

   Description:
   ============

   Calculate RMS deviation between each conformation and the reference 
   set in a CSearch (CHARMM-free CONGEN) CG file.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  01.06.94 Original
   V1.1  27.09.94 Added Usage message

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"

#include "cgutil.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL HasDummyAtoms(COOR *conf, HEADER header);
REAL CalcRMS(COOR *RefConf, COOR *conf, int ncoor);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, int *type, char *RefPDB, 
                  char *CGFile);
COOR *GetRequiredAtoms(COOR *conf, HEADER header, int type, PDB *pdb, 
                       int *natoms);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for counting conformations

   01.06.94 Original    By: ACRM
   27.09.94 Added -h check
*/
int main(int argc, char **argv)
{
   HEADER header;
   COOR   *conf;
   REAL   energy;
   FILE   *in = stdin;
   int    nconf;

   if(argc > 1)
   {
      if(!strncmp(argv[1],"-h",2) || !strncmp(argv[1],"-H",2))
      {
         Usage();
         return(0);
      }
      
      if((in = fopen(argv[1],"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",argv[1]);
         return(1);
      }
   }

   if(!ReadHeader(in, &header))
   {
      fprintf(stderr,"File has bad header\n");
      return(1);
   }
   else
   {
      conf = GetConf(in, header, &energy);
      if(conf == NULL || conf == (COOR *)(-1))
      {
         fprintf(stderr,"Error in reading reference conformation\n");
         return(1);
      }

      nconf = 0;
      while((conf = GetConf(in, header, &energy))!=NULL)
      {
         if(conf == (COOR *)(-1))
         {
            fprintf(stderr,"Error in conformation file\n");
            fprintf(stderr,"%d conformations read OK\n",nconf);
            
            return(1);
         }
         free(conf);
         nconf++;
      }

      printf("%d\n",nconf);
   }
   
   return(0);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   27.09.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nCountCG V1.1 (C) 1994 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: countcg [file.cg]\n");
   fprintf(stderr,"(Input is read from stdin if no file specified.)\n\n");
   fprintf(stderr,"Counts the conformations in a Charmm-free CONGEN CG \
file, ignoring the\n");
   fprintf(stderr,"reference set.\n\n");
}
