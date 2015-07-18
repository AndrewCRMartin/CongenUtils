/*************************************************************************

   Program:    patchcg
   File:       patchcg.c
   
   Version:    V1.1
   Date:       26.05.94
   Function:   Patch the reference set of a CG file with the true coords
   
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
   V1.0  26.05.94 Original
   V1.1  10.11.94 Added warning message about atom order

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
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
BOOL WritePDBReferenceSet(FILE *fp, PDB *pdb, HEADER header);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for paching CSearch .CG files with reference coords

   26.05.94 Original    By: ACRM
   10.11.94 Added warning message.
*/
int main(int argc, char **argv)
{
   HEADER header;
   FILE   *cgin_fp,
          *cgout_fp,
          *pdb_fp;
   COOR   *conf;
   PDB    *pdb;
   REAL   energy;
   int    nconf,
          natoms;
      
   if(argc != 4)
   {
      Usage();
      exit(1);
   }

   /* Open files                                                        */
   if((cgin_fp = fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input CG file; %s\n",argv[1]);
      return(1);
   }
   if((pdb_fp = fopen(argv[2],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input PDB file; %s\n",argv[2]);
      return(1);
   }
   if((cgout_fp = fopen(argv[3],"w"))==NULL)
   {
      fprintf(stderr,"Unable to open output CG file; %s\n",argv[3]);
      return(1);
   }

   /* Read the header from the input CG file                            */
   if(!ReadHeader(cgin_fp, &header))
   {
      fprintf(stderr,"CG file has bad header\n");
      return(1);
   }

   /* Read the PDB file                                                 */
   if((pdb = ReadPDB(pdb_fp, &natoms))==NULL)
   {
      fprintf(stderr,"Unable to read PDB file\n");
      return(1);
   }

   /* Write header to the output CG file                                */
   WriteHeader(cgout_fp, header);

   /* Read the reference conf from the CG file                          */
   if((conf = GetConf(cgin_fp, header, &energy))==NULL)
   {
      fprintf(stderr,"Input CG file has no data!\n");
      return(1);
   }
   free(conf);
   
   /* Get the reference atom coords from the PDB and write them out     */
   if(!WritePDBReferenceSet(cgout_fp, pdb, header))
   {
      fprintf(stderr,"Unable to find all reference atoms in the PDB \
file\n");
      return(1);
   }
   
   /* Now copy remaining conformations                                  */
   nconf = 0;
   while((conf = GetConf(cgin_fp, header, &energy))!=NULL)
   {
      nconf++;
      if(conf==((COOR *)(-1)))
      {
         fprintf(stderr,"Error in input conformation file at \
conf. %d\n",nconf);
         return(1);
      }
      PutConf(cgout_fp,header,conf,energy);
      free(conf);
   }

   /* Close files, free memory and exit                                 */
   fclose(cgin_fp);
   fclose(cgout_fp);
   fclose(pdb_fp);
   KillHeader(&header);

   /* Print warning message                                             */
   fprintf(stderr,"N.B. Results will be wrong if the PDB file does not \
have the correct\n");
   fprintf(stderr,"atom ordering and hydrogens (i.e. N, H, CA, s/c, C, \
O)\n");
   
   
   return(0);
}

/************************************************************************/
/*>BOOL WritePDBReferenceSet(FILE *fp, PDB *pdb, HEADER header)
   ------------------------------------------------------------
   Input:   FILE    *fp        Output CG file pointer
            PDB     *pdb       PDB linked list from which to get
                               reference atoms
            HEADER  header     CG header structure
   Returns: BOOL               Success (all atoms found)

   Searches the PDB linked list for the atom numbers specified in the
   header structure and writes them to the output file. Finally writes
   an enrgy value of 0.0. Returns TRUE if all atoms found; else FALSE.

   26.05.94 Original    By: ACRM
*/
BOOL WritePDBReferenceSet(FILE *fp, PDB *pdb, HEADER header)
{
   int  i;
   PDB  *p, *pos;
   BOOL found;
   
   p = pos = pdb;
   
   for(i=0; i<header.ncons; i++)
   {
      found = FALSE;
      
      /* Search on from pos to the end of the linked list               */
      for(p=pos; p!=NULL; NEXT(p))
      {
         if(p->atnum == header.ConsAtoms[i])
         {
            found = TRUE;
            break;
         }
      }
      
      /* If not found, search from beginning of list to pos             */
      if(!found)
      {
         for(p=pdb; p!=pos && p!=NULL; NEXT(p))
         {
            if(p->atnum == header.ConsAtoms[i])
            {
               found = TRUE;
               break;
            }
         }
      }
      
      /* If still not found, we have an error                           */
      if(!found)
         return(FALSE);
      
      /* Found, so write coords, update pos and continue                */
      fprintf(fp,"%8.3f %8.3f %8.3f\n", p->x, p->y, p->z);
      pos = p;
   }

   fprintf(fp,"%15.5f\n",(double)0.0);
   return(TRUE);
}



/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a Usage message

   26.05.94 Original    By: ACRM
   10.11.94 Added warnings about atom order
*/
void Usage(void)
{
   fprintf(stderr,"\npatchcg  V1.1\n"); 
   fprintf(stderr,"=============\n\n");
   fprintf(stderr,"(c) 1994 Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in\nso doing.\n\n");
   fprintf(stderr,"Adds a true reference coordinate set to a CSearch \
(Charmm-free CONGEN) \nconformation file. Coords are extracted from \
the PDB file (with \nhydrogens) used by CSearch\n\n");
   fprintf(stderr,"Usage: patchcg <input.cg> <reference.pdb> \
<output.cg>\n\n");
   fprintf(stderr,"N.B. The PDB file *must* have hydrogens and correct \
atom ordering\n");
   fprintf(stderr,"(i.e. N, H, CA, s/c, C, O)\n");
   fprintf(stderr,"No warning is given if this is not the case!!\n");
   fprintf(stderr,"You can use pdborder -c -i to get the correct \
order.\n\n");
}
