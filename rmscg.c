/*************************************************************************

   Program:    rmscg
   File:       rmscg.c
   
   Version:    V1.1
   Date:       27.09.94
   Function:   Calculate RMS deviation between each conformation and
               the reference set.
   
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
   V1.0  26.05.94 Original
   V1.1  27.09.94 Stopped crash on rmscg -h

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
#define TYPE_ALL    0     /* RMS atom types to include                  */
#define TYPE_NOH    1
#define TYPE_NCACO  2
#define TYPE_CA     3

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
/*>COOR *GetRequiredAtoms(COOR *conf, HEADER header, int type, PDB *pdb, 
                          int *natoms)
   ---------------------------------------------------------------------
   Input:   COOR    *conf       Coordinate array
            HEADER  header      CG header structure
            int     type        Atom selection type (TYPE_NOH, etc)
            PDB     *pdb        PDB linked list from which to identify
                                atom types
   I/O:     int     *natoms     Number of atoms before & after selection
   Returns: COOR    *           Coordinate array after selecting atoms

   Takes a coordinate array (which has been created by malloc()), a
   CG header and a PDB linked list and creates a new coordinate array
   containing only atoms of the required types. The original array
   is freed.

   26.05.94 Original    By: ACRM
   27.05.94 Added H case
*/
COOR *GetRequiredAtoms(COOR *conf, HEADER header, int type, PDB *pdb, 
                       int *natoms)
{
   COOR       *retconf;
   static int *atomlist = NULL,
              NumAtoms;
   int        i, j;
   PDB        *p;
   
   /* If this is the first call, create a list of atoms we need         */
   if(atomlist == NULL)
   {
      /* Allocate memory for the list                                   */
      if((atomlist = malloc((*natoms) * sizeof(int)))==NULL)
      {
         *natoms = (-1);
         return(conf);
      }

      NumAtoms = 0;

      /* Go through the complete list                                   */
      for(i=0; i<header.ncons; i++)
      {
         /* Run through the PDB linked list                             */
         for(p=pdb; p!=NULL; NEXT(p))
         {
            /* If we get the atom number                                */
            if(p->atnum == header.ConsAtoms[i])
            {
               /* Copy in the atom's number if it's a required atom     */
               switch(type)
               {
               case TYPE_NCACO:
                  if(!strncmp(p->atnam,"C   ",4) ||
                     !strncmp(p->atnam,"N   ",4) ||
                     !strncmp(p->atnam,"O   ",4))
                  {
                     atomlist[NumAtoms++] = p->atnum;
#ifdef DEBUG
                     WritePDBRecord(stderr,p);
#endif
                  }
                  /* Fall through......                                 */
               case TYPE_CA:
                  if(!strncmp(p->atnam,"CA  ",4))
                  {
                     atomlist[NumAtoms++] = p->atnum;
#ifdef DEBUG
                     WritePDBRecord(stderr,p);
#endif
                  }
                  break;
               case TYPE_NOH:
                  if(p->atnam[0] != 'H')
                  {
                     atomlist[NumAtoms++] = p->atnum;
#ifdef DEBUG
                     WritePDBRecord(stderr,p);
#endif
                  }
                  break;
               default:
                  atomlist[NumAtoms++] = p->atnum;
#ifdef DEBUG
                  WritePDBRecord(stderr,p);
#endif
               }
            }
         }
      }
   }

   *natoms = NumAtoms;
   
   /* Allocate space for the output coordinate array                    */
   if((retconf = (COOR *)malloc(NumAtoms * sizeof(COOR)))==NULL)
   {
      *natoms = (-1);
      return(conf);
   }

   /* Copy in coords of required atoms                                  */
   for(i=0, j=0; i<header.ncons; i++)
   {
      if(header.ConsAtoms[i] == atomlist[j])
      {
         retconf[j].x = conf[i].x;
         retconf[j].y = conf[i].y;
         retconf[j].z = conf[i].z;
         j++;
      }
   }

   if(j != NumAtoms)
   {
      fprintf(stderr,"Programmer error in GetRequiredAtoms()\n");
      exit(1);
   }
   
   free(conf);
   return(retconf);
}


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for calculating RMS of each conf in a CSearch CG file

   24.05.94 Original    By: ACRM
   26.05.94 Added call to ParseCmdLine() & GetRequiredAtoms() and PDB
            reading stuff.
*/
int main(int argc, char **argv)
{
   HEADER header;
   COOR   *RefConf,
          *conf;
   FILE   *in;
   REAL   energy,
          rms;
   int    nconf,
          natoms,
          type;
   char   RefPDB[160],
          CGFile[160];
   PDB    *pdb;

   if(!ParseCmdLine(argc, argv, &type, RefPDB, CGFile))
   {
      Usage();
      return(1);
   }

   if((in = fopen(CGFile,"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file: %s\n",CGFile);
      return(1);
   }

   if(!ReadHeader(in, &header))
   {
      fprintf(stderr,"File has bad header\n");
      return(1);
   }
   else
   {
      if(type != TYPE_ALL)
      {
         FILE *fp;
         
         if((fp=fopen(RefPDB,"r"))==NULL)
         {
            fprintf(stderr,"Unable to open reference PDB file: %s\n",
                    RefPDB);
            return(1);
         }

         if((pdb = ReadPDB(fp, &natoms))==NULL)
         {
            fprintf(stderr,"No atoms read from reference PDB file\n");
            return(1);
         }
      }
      
      RefConf = GetConf(in, header, &energy);
      if(RefConf == NULL || RefConf == (COOR *)(-1))
      {
         fprintf(stderr,"Error in reading reference conformation\n");
         return(1);
      }

      if(type != TYPE_ALL)
      {
         RefConf = GetRequiredAtoms(RefConf, header, type, pdb, &natoms);
         if(natoms == (-1))
         {
            fprintf(stderr,"No memory to get required atoms\n");
            return(1);
         }
#ifdef DEBUG
         fprintf(stderr,"Atoms remaining: %d\n",natoms);
#endif
      }

      if(HasDummyAtoms(RefConf,header))
      {
         fprintf(stderr,"Reference set contains dummy atoms.\n");
         return(0);
      }
      
      printf("Reference set energy: %f\n",energy);

      nconf = 1;
      while((conf = GetConf(in, header, &energy))!=NULL)
      {
         if(conf == (COOR *)(-1))
         {
            fprintf(stderr,"Error in conformation file\n");
            return(1);
         }

         natoms = header.ncons;

         if(type != TYPE_ALL)
         {
            conf = GetRequiredAtoms(conf, header, type, pdb, &natoms);
            if(natoms == (-1))
            {
               fprintf(stderr,"No memory to get required atoms\n");
               return(1);
            }
         }
         
         rms = CalcRMS(RefConf, conf, natoms);
         printf("Conformation %5d: RMS: %8.3f Energy: %10.5f\n",
                nconf++, rms, energy);
         free(conf);
      }
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL HasDummyAtoms(COOR *conf, HEADER header)
   ---------------------------------------------
   Input:   COOR   *conf     Coordinate array
            HEADER header    CG header structure
   Returns: BOOL             Has dummy atoms?

   Checks whether a conformation stored in a coordinate array has any
   dummy atoms.

   24.05.94 Original    By: ACRM
*/
BOOL HasDummyAtoms(COOR *conf, HEADER header)
{
   int i;
   
   for(i=0; i<header.ncons; i++)
   {
      if(conf[i].x > (REAL)9998.0 &&
         conf[i].y > (REAL)9998.0 &&
         conf[i].z > (REAL)9998.0)
         return(TRUE);
   }
   
   return(FALSE);
}

/************************************************************************/
/*>REAL CalcRMS(COOR *RefConf, COOR *conf, int ncoor)
   --------------------------------------------------
   Input:   COOR     *RefCoor     Reference coordinate array
            COOR     *conf        Comparison coordinate array
            int      ncoor        Number of coordinates to compare

   Calculates the RMS of two COOR arrays.

   24.05.94 Original    By: ACRM
*/
REAL CalcRMS(COOR *RefConf, COOR *conf, int ncoor)
{
   int  i;
   REAL rms;
   
   for(i=0, rms = (REAL)0.0; i<ncoor; i++)
      rms += DISTSQ((&(RefConf[i])), (&(conf[i])));

   rms /= ncoor;

   rms = (REAL)sqrt((double)rms);
   
   return(rms);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, int *type, char *RefPDB, 
                     char *CGFile)
   -----------------------------------------------------------------
   Input:   int     argc     Argument count
            char    **argv   Argument array
   Output:  int     *type    RMS type (TYPE_ALL, TYPE_CA, TYPE_NCACO)
            char    *RefPDB  Reference PDB (if TYPE_CA or TYPE_NCACO)
            char    *CGFile  CG file name
   Returns: BOOL             Success?

   Parses the command line.

   26.05.94 Original    By: ACRM
   27.05.94 Added H case
   27.09.94 Added check on arg count in while loop
*/
BOOL ParseCmdLine(int argc, char **argv, int *type, char *RefPDB, 
                  char *CGFile)
{
   argc--;
   argv++;

   *type   = TYPE_ALL;
   *RefPDB = '\0';

   while(argc && argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 'h':
      case 'H':
         *type = TYPE_NOH;
         argv++;
         argc--;
         if(argc) strcpy(RefPDB,argv[0]);
         break;
      case 'c':
      case 'C':
         *type = TYPE_CA;
         argv++;
         argc--;
         if(argc) strcpy(RefPDB,argv[0]);
         break;
      case 'b':
      case 'B':
         *type = TYPE_NCACO;
         argv++;
         argc--;
         if(argc) strcpy(RefPDB,argv[0]);
         break;
      default:
         return(FALSE);
      }
      argv++;
      argc--;
   }

   if(argc != 1)
      return(FALSE);
   
   strcpy(CGFile, argv[0]);
   return(TRUE);
}
      
/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   24.05.94 Original    By: ACRM
   27.05.94 Added H info
*/
void Usage(void)
{
   fprintf(stderr,"\nrmscg  V1.1\n"); 
   fprintf(stderr,"===========\n\n");
   fprintf(stderr,"(c) 1994 Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in\nso doing.\n\n");
   fprintf(stderr,"Calculates the RMS of each conformation in a CSearch \
(Charmm-free CONGEN) \nconformation files compared with the \
reference set.\n\n");
   fprintf(stderr,"Usage: rmscg [-h <ref.pdb>] [-c <ref.pdb>] \
[-b <ref.pdb>] <input.cg>\n");
   fprintf(stderr,"       -h Calculate All atom, but no H RMS\n");
   fprintf(stderr,"       -c Calculate CA RMS\n");
   fprintf(stderr,"       -b Calculate N,CA,C,O RMS\n");
}

