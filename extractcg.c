/*************************************************************************

   Program:    extractcg
   File:       extractcg.c
   
   Version:    V1.1
   Date:       12.10.94
   Function:   Patch a coordinate file with a conformation from a CG file
   
   Copyright:  (c) Dr. Andrew C. R. Martin, 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling,
               University College London,
   Phone:      HOME: +44 (0)372 275775
   EMail:      JANET:    martin@uk.ac.ucl.bioc.bsm
               INTERNET: martin@bsm.bioc.ucl.ac.uk
               
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
   V1.0  21.09.94 Original    By: ACRM
   V1.1  12.10.94 Fixed bug in command line parsing

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "cgutil.h"

/************************************************************************/
/* Defines
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void DoProcessing(FILE *out, FILE *cgfp, PDB *pdb, int ConfNum);
void PatchPDBWithConf(PDB *pdb, COOR *conf, HEADER header);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, int *ConfNum, char *CGFile, 
                  char *InPDBFile, char *OutPDBFile);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for paching CSearch .CG files with reference coords

   21.09.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   int    ConfNum = 1,
          natoms;
   FILE   *in     = stdin,
          *out    = stdout,
          *cgfp   = NULL;
   char   InPDBFile[MAXBUFF],
          OutPDBFile[MAXBUFF],
          CGFile[MAXBUFF];
   PDB    *pdb;
   
      
   if(ParseCmdLine(argc, argv, &ConfNum, CGFile, InPDBFile, OutPDBFile))
   {
      if(OpenStdFiles(InPDBFile, OutPDBFile, &in, &out))
      {
         /* Open the CG file                                            */
         if((cgfp = fopen(CGFile,"r"))!=NULL)
         {
            /* Read in the PDB file                                     */
            if((pdb = ReadPDB(in, &natoms))!=NULL)
            {
               DoProcessing(out, cgfp, pdb, ConfNum);
               FREELIST(pdb, PDB);
            }
            else
            {
               fprintf(stderr,"Unable to read PDB file\n");
               return(1);
            }

            fclose(cgfp);
         }
         else
         {
            fprintf(stderr,"Unable to open conformation file: %s\n",
                    CGFile);
            return(1);
         }
      }
   }
   else
   {
      Usage();
      return(1);
   }

   return(0);
}

/************************************************************************/
/*>void DoProcessing(FILE *out, FILE *cgfp, PDB *pdb, int ConfNum)
   ---------------------------------------------------------------
   Does the actual processing work and writes a new PDB file.
   The required conformation is extracted from the conformation file,
   its coordinates are patched into the PDB linked list and the new
   file is written.

   21.09.94 Original    By: ACRM
*/
void DoProcessing(FILE *out, FILE *cgfp, PDB *pdb, int ConfNum)
{
   HEADER header;
   int    i;
   COOR   *conf = NULL;
   REAL   energy;

   /* Read the header from the input CG file                            */
   if(!ReadHeader(cgfp, &header))
   {
      fprintf(stderr,"Conformation file has bad header\n");
      return;
   }

   /* Read along to the appropriate conformation                        */
   for(i=0; i<=ConfNum; i++)
   {
      /* If we've already read a conformation (i.e. not the first time
         around the loop), free the memory allocated.
      */
      if(conf != NULL)
         free(conf);
      
      /* Read the reference conf from the CG file                       */
      if((conf = GetConf(cgfp, header, &energy))==NULL)
      {
         fprintf(stderr,"Conformation number %d not found\n",ConfNum);
         return;
      }
   }
   
   /* Patch the coordinates into the PDB linked list                    */
   PatchPDBWithConf(pdb, conf, header);

   /* Write out the patched PDB file                                    */
   WritePDB(out, pdb);

   /* Free memory and return                                            */
   if(conf != NULL) free(conf);
   KillHeader(&header);
}

/************************************************************************/
/*>void PatchPDBWithConf(PDB *pdb, COOR *conf, HEADER header)
   ----------------------------------------------------------
   Patches the coordinates from conf into the PDB linked list

   21.09.94 Original    By: ACRM
*/
void PatchPDBWithConf(PDB *pdb, COOR *conf, HEADER header)
{
   int i,
       AtomNum;
   PDB *p;
   
   /* Initialise to point to the first atom in the linked list          */
   p       = pdb;
   AtomNum = 1;

   for(i=0; i<header.ncons; i++)
   {
      /* If the current atom number is beyond the one we're looking for,
         go back to the beginning.
      */
      if(AtomNum > header.ConsAtoms[i])
      {
         p       = pdb;
         AtomNum = 1;
      }
      else if(AtomNum < header.ConsAtoms[i])
      {
         /* Current atom is before the one we are looking for, so step on
            to the appropriate atom
         */
         while(AtomNum < header.ConsAtoms[i])
         {
            if(p!=NULL) NEXT(p);
            AtomNum++;
         }
      }

      if(p!=NULL)
      {
         /* We should now be at the correct atom, so patch it           */
         p->x = conf[i].x;
         p->y = conf[i].y;
         p->z = conf[i].z;
      }
      else
      {
         /* The atom number went off the end of the list, so reset      */
         p       = pdb;
         AtomNum = 1;
      }
   }
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a Usage message

   21.09.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nextractcg  V1.1\n"); 
   fprintf(stderr,"===============\n\n");
   fprintf(stderr,"(c) 1994 Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in\nso doing.\n\n");
   fprintf(stderr,"Extracts a conformation from a CSearch (Charmm-free \
CONGEN) conformation\n");
   fprintf(stderr,"file. The PDB file must have NTER and CTER residues \
and hydrogens as\n");
   fprintf(stderr,"required by CSearch\n\n");

   fprintf(stderr,"Usage: extractcg [-<n>] <input.cg> [<input.pdb>] \
[<output.cg>]\n\n");

   fprintf(stderr,"If either input or output file is not specified, \
standard input/output\n");
   fprintf(stderr,"are used. If the conformation number is not specified, \
the first\n");
   fprintf(stderr,"conformation will be extracted.\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, int *ConfNum, char *CGFile, 
                     char *InPDBFile, char *OutPDBFile)
   --------------------------------------------------------------------
   Parse the command line. Requires the CGfile to be specified; other
   parameters are optional.

   21.09.94 Original    By: ACRM
   12.10.94 Changed !isalpha() to !isdigit
*/
BOOL ParseCmdLine(int argc, char **argv, int *ConfNum, char *CGFile, 
                  char *InPDBFile, char *OutPDBFile)
{
   argc--;
   argv++;

   InPDBFile[0] = OutPDBFile[0] = CGFile[0] = '\0';

   if(argc == 0)
      return(FALSE);

   while(argv[0][0] == '-')
   {
      if(!isdigit(argv[0][1]))
      {
         return(FALSE);
      }
      if(sscanf(argv[0]+1,"%d",ConfNum) == 0)
         return(FALSE);
      
      argc--;
      argv++;
   }
   
   /* There must be at least one argument                               */
   if(argc < 1 || argc > 3)
      return(FALSE);
   strcpy(CGFile, argv[0]);
   argc--;
   argv++;
   
   /* If there's another argument, copy it                              */
   if(argc)
   {
      strcpy(InPDBFile, argv[0]);
      argc--;
      argv++;
   }
   
   /* If there's another argument, copy it                              */
   if(argc)
   {
      strcpy(OutPDBFile, argv[0]);
      argc--;
      argv++;
   }
   
   return(TRUE);
}

