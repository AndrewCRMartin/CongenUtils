/*************************************************************************

   Program:    spancg
   File:       spancg.c
   
   Version:    V1.0
   Date:       15.09.94
   Function:   Scan a CG file and find the greatest movement of any atom
               from the reference coordinate position
   
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
void Usage(void);

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message.

   15.09.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nspancg  V1.0\n"); 
   fprintf(stderr,"============\n\n");
   fprintf(stderr,"(c) 1994 Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in\nso doing.\n\n");
   fprintf(stderr,"Finds the greatest atom movement from the reference \
set coordinates\n");

   fprintf(stderr,"in a CSearch (Charmm-free CONGEN) conformation \
file.\n\n");
   fprintf(stderr,"Usage: spancg <input.cg>\n\n");
}

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for CSearch span measurement

   15.09.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   int    i;
   HEADER header;
   FILE   *in;
   COOR   *coor    = NULL,
          *refcoor = NULL;
   REAL   energy,
          dist,
          maxdist = (REAL)0.0;
   
   /* Check the command line                                            */
   if(argc < 2)
   {
      Usage();
      return(1);
   }

   /* Open the input file                                               */
   if((in=fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file: %s\n",argv[1]);
      return(1);
   }

   /* Read input file header                                            */
   if(!ReadHeader(in,&header))
   {
      fprintf(stderr,"File %s has a bad header\n",argv[1]);
      return(1);
   }

   while((coor=GetConf(in,header,&energy))!=NULL)
   {
      /* Check for error return                                         */
      if(coor == (COOR *)(-1))
      {
         fprintf(stderr,"No memory to read coordinates\n");
         return(1);
      }

      if(refcoor == NULL)   /* Just keep the reference coordinates      */
      {
         refcoor = coor;
         for(i=0; i<header.ncons; i++)
         {
            if(coor[i].x > (REAL)9998.0 &&
               coor[i].y > (REAL)9998.0 &&
               coor[i].z > (REAL)9998.0)
            {
               fprintf(stderr,"Reference set must not contain dummy \
coordinates. Use patchcg.\n");
               return(1);
            }
         }
      }
      else                  /* Find largest movement                    */
      {
         for(i=0; i<header.ncons; i++)
         {
            dist = DIST(&(refcoor[i]), &(coor[i]));
            if(dist > maxdist)
               maxdist = dist;
         }
         free(coor);
      }
   }

   if(refcoor != NULL) free(refcoor);

   printf("Maximum atom movement was %f Angstroms\n",maxdist);

   return(0);
}

