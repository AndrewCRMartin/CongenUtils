/*************************************************************************

   Program:    mergecg
   File:       mergecg.c
   
   Version:    V1.0
   Date:       24.05.94
   Function:   Merge 2 CSearch .cg files
   
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
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================


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
int main(int argc, char **argv);
int CopyConfs(FILE *in, FILE *out, HEADER header);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for CSearch CG merging

   24.05.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   int    nconf,
          i,
          TotalConf = 0;
   HEADER header1,
          header2;
   FILE   *in1, *in2, *out;
   
   /* Check the command line                                            */
   if(argc < 3)
   {
      Usage();
      return(1);
   }

   /* Open the first input file and the output file                     */
   if((in1=fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file: %s\n",argv[1]);
      return(1);
   }
   if((out=fopen(argv[argc-1],"w"))==NULL)
   {
      fprintf(stderr,"Unable to open output file: %s\n",argv[argc-1]);
      return(1);
   }

   /* Read first input file header                                      */
   if(!ReadHeader(in1,&header1))
   {
      fprintf(stderr,"File %s has a bad header\n",argv[1]);
      return(1);
   }
   D("Read Header\n");

   /* Write the header to the output file                               */
   WriteHeader(out,header1);
   D("Written Header\n");
   
   /* Copy the conformations over to the output file and close the input*/
   nconf = CopyConfs(in1,out,header1) - 1;
   fclose(in1);
   D("Copied first file confs\n");

   /* Check for allocation failure                                      */
   if(nconf == (-1))
   {
      fprintf(stderr,"No memory for conformation\n");
      return(1);
   }

   /* Print number of conformations                                     */
   printf("%s: %d conformations.\n",argv[1],nconf);
   TotalConf = nconf;
   
   /* For each remaining input file                                     */
   for(i=2; i<argc-1; i++)
   {
      /* Open the file                                                  */
      if((in2=fopen(argv[i],"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",argv[i]);
         return(1);
      }

      /* Read the header                                                */
      if(!ReadHeader(in2,&header2))
      {
         fprintf(stderr,"File %s has a bad header\n",argv[i]);
         return(1);
      }
      
      /* Check it matches the first header                              */
      if(HeadersMatch(header1,header2))
      {
         COOR *coor;
         REAL energy;
         
         /* Get and skip the reference conformation                     */
         coor=GetConf(in2,header1,&energy);
         if(coor != NULL) free(coor);

         /* Copy the conformations over to the output file and close the 
            input
         */
         nconf = CopyConfs(in2,out,header1);
         fclose(in2);

         /* Check for allocation failure                                */
         if(nconf == (-1))
         {
            fprintf(stderr,"No memory for conformation\n");
            return(1);
         }

         KillHeader(&header2);

         /* Print number of conformations                               */
         printf("%s: %d conformations.\n",argv[i],nconf);
         TotalConf += nconf;
      }
      else   /* Header mismatch                                         */
      {
         fprintf(stderr,"Header on file %s does not match (skipped)\n", 
                 argv[i]);
      }
   }

   printf("Total conformations: %d\n",TotalConf);

   return(0);
}

/************************************************************************/
/*>int CopyConfs(FILE *in, FILE *out, HEADER header)
   -------------------------------------------------
   Input:   FILE   *in      Input CG file pointer
            FILE   *out     Output CG file pointer
            HEADER header   CG file header structure

   Copies call remaining onformations from one CG file to another. 
   Assumes that the headers have been read from each and that one header 
   is stored in the header structure.

   24.05.94 Original    By: ACRM
*/
int CopyConfs(FILE *in, FILE *out, HEADER header)
{
   COOR *coor;
   REAL energy;
   
   int nconf=0;

   while((coor=GetConf(in,header,&energy))!=NULL)
   {
      /* Check for error return                                         */
      if(coor == (COOR *)(-1))
         return(-1);

      D("Got conf.\n");

      /* Write the conformation                                         */
      PutConf(out,header,coor,energy);
      
      nconf++;
   }

   free(coor);
   return(nconf);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message.

   24.05.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nmergecg  V1.0\n"); 
   fprintf(stderr,"=============\n\n");
   fprintf(stderr,"(c) 1994 Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in\nso doing.\n\n");
   fprintf(stderr,"Merges any number of CSearch (Charmm-free CONGEN) \
conformation files.\n\n");
   fprintf(stderr,"Usage: mergecg <input.cg> [<input2.cg> ...] \
<output.cg>\n\n");
}
