/*************************************************************************

   Program:    
   File:       cgutil.c
   
   Version:    V1.0
   Date:       24.05.94
   Function:   General purpose utility routines for handling CSearch 
               CG files
   
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
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "cgutil.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/*>BOOL ReadHeader(FILE *fp, HEADER *header)
   -----------------------------------------
   Input:   FILE   *fp      CG file pointer
   Output:  HEADER *header  Completed header structure

   Reads a header from a CG file. memory is allocated within the header 
   structure to store the constructued atoms list

   24.05.94 Original    By: ACRM
*/
BOOL ReadHeader(FILE *fp, HEADER *header)
{
   char buffer[160],
        word[16],
        *bp;
   int  ncons;
   BOOL GotValues = FALSE;

   /* First line contains nres, natoms, ncons                           */
   fgets(buffer,159,fp);
   TERMINATE(buffer);
   
   if((sscanf(buffer,"%d %d %d",
              &(header->nres),&(header->natoms),&(header->ncons))) != 3)
      return(FALSE);

   /* Allocate memory for constructed atoms array                       */
   if((header->ConsAtoms = (int *)malloc(header->ncons * sizeof(int)))
      ==NULL)
      return(FALSE);

   /* Read the constructed atom data                                    */
   for(ncons = 0, GotValues = FALSE; ncons < header->ncons; )
   {
      if(!GotValues)
      {
         fgets(buffer,159,fp);
         TERMINATE(buffer);
         GotValues = TRUE;
         bp = buffer;
      }

      bp = GetWord(bp, word);
      if(strlen(word))
      {
         if(!sscanf(word,"%d",&(header->ConsAtoms[ncons])))
            return(FALSE);
         ncons++;
      }
      else
      {
         GotValues = FALSE;
      }
   }
   
   return(TRUE);
}
   
/************************************************************************/
/*>void WriteHeader(FILE *fp, HEADER header)
   -----------------------------------------
   Input:   FILE    *fp     Output CG file pointer
            HEADER  header  CG file header structure

   Writes a CG header

   24.05.94 Original    By: ACRM
*/
void WriteHeader(FILE *fp, HEADER header)
{
   int i;
   
   fprintf(fp,"%5d %5d %5d", header.nres, header.natoms, header.ncons);
   for(i=0; i<header.ncons; i++)
   {
      if(!(i%13))
         fprintf(fp,"\n");

      fprintf(fp,"%5d ",header.ConsAtoms[i]);
   }
   fprintf(fp,"\n");
}

/************************************************************************/
/*>void KillHeader(HEADER *header)
   -------------------------------
   Output:  HEADER   *header     Pointer to a header structure

   Clears the values in a header structure and frees memory allocated
   for the constructed atoms list

   24.05.94 Original    By: ACRM
*/
void KillHeader(HEADER *header)
{
   header->natoms = 0;
   header->nres   = 0;
   header->ncons  = 0;
   if(header->ConsAtoms != NULL)
      free(header->ConsAtoms);
   header->ConsAtoms = NULL;
}
   
/************************************************************************/
/*>COOR *GetConf(FILE *fp, HEADER header, REAL *energy)
   ----------------------------------------------------
   Input:   FILE    *fp       Input CG file pointer
            HEADER  header    Header structure from CG file
   Output:  REAL    *energy   Energy for this conformation
   Returns: COOR    *         Allocated array of coordinates for this
                              conformation
                              NULL if no more confs
                              -1   if error

   Reads the next conformation from a CG file. Allocates memory for
   the coordinates.

   24.05.94 Original    By: ACRM
*/
COOR *GetConf(FILE *fp, HEADER header, REAL *energy)
{
   int  i;
   COOR *coor;
   BOOL error = FALSE;
   char buffer[160];
   
   if((coor = (COOR *)malloc(header.ncons * sizeof(COOR)))==NULL)
      return((COOR *)(-1));
   
   for(i=0; i<header.ncons; i++)
   {
      if((fgets(buffer,159,fp)))
      {
         TERMINATE(buffer);
         
         if((sscanf(buffer,"%lf %lf %lf",
                    &(coor[i].x), &(coor[i].y), &(coor[i].z)))!=3)
         {
            D("Error in coordinate line: ");
            D(buffer);
            D("\n");
            
            error = TRUE;
            break;
         }
      }
      else
      {
         if(i!=0)
         {
            D("Missing coordinate line\n");
            
            error = TRUE;
            break;
         }
         else
         {
            free(coor);
            return(NULL);
         }
      }
   }

   if(!error)
   {
      if((fgets(buffer,159,fp)))
      {
         if((sscanf(buffer,"%lf",energy))!=1)
            error = TRUE;
      }
      else
      {
         error = TRUE;
      }
   }

   if(error)
   {
      free(coor);
      return((COOR *)(-1));
   }
   
   return(coor);
}

/************************************************************************/
/*>void PutConf(FILE *out, HEADER header, COOR *coor, REAL energy)
   ---------------------------------------------------------------
   Input:   FILE    *out     Output CG file pointer
            HEADER  header   CG header structure
            COOR    *coor    Array of coordinates for conformation
            REAL    energy   Energy for this conformation

   Writes a conformation to an output CG file

   24.05.94 Original    By: ACRM
*/
void PutConf(FILE *out, HEADER header, COOR *coor, REAL energy)
{
   int i;
   
   for(i=0; i<header.ncons; i++)
      fprintf(out,"%8.3f %8.3f %8.3f\n",coor[i].x,coor[i].y,coor[i].z);

   fprintf(out,"%15.5f\n",energy);
}

/************************************************************************/
/*>BOOL HeadersMatch(HEADER header1, HEADER header2)
   -------------------------------------------------
   Input:   HEADER    header1     First header structure
            HEADER    header2     Second header structure
   Returns: BOOL                  Match?

   Checks whether 2 CG headers match.

   24.05.94 Original    By: ACRM
*/
BOOL HeadersMatch(HEADER header1, HEADER header2)
{
   int i;
   
   if(header1.nres   != header2.nres)   return(FALSE);
   if(header1.natoms != header2.natoms) return(FALSE);
   if(header1.ncons  != header2.ncons)  return(FALSE);

   for(i=0; i<header1.ncons; i++)
   {
      if(header1.ConsAtoms[i] != header2.ConsAtoms[i]) return(FALSE);
   }
   
   return(TRUE);
}


