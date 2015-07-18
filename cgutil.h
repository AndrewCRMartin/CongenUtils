/*************************************************************************

   Program:    
   File:       cgutil.h
   
   Version:    V1.0
   Date:       24.05.94
   Function:   Header file for CG utility routines
   
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
#ifndef _CGUTIL_H
#define _CGUTIL_H

/* Includes
*/
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"

/************************************************************************/
/* Defines
*/
typedef struct
{
   int *ConsAtoms,
       nres,
       natoms,
       ncons;
}  HEADER;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ReadHeader(FILE *fp, HEADER *header);
void WriteHeader(FILE *fp, HEADER header);
void KillHeader(HEADER *header);
COOR *GetConf(FILE *fp, HEADER header, REAL *energy);
void PutConf(FILE *out, HEADER header, COOR *coor, REAL energy);
BOOL HeadersMatch(HEADER header1, HEADER header2);



#endif /* _CGUTIL_H                                                     */
