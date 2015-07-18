/*************************************************************************

   Program:    ClusterCG
   File:       clustercg.c
   
   Version:    V1.1
   Date:       01.02.95
   Function:   Cluster conformations from CG file
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994-5
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

   clustercg [-n target] [-a angle] [-m maxang] in.cg ref.pdb out.cg

**************************************************************************

   Revision History:
   =================
   V1.0  30.11.94 Original    By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/MathType.h"
#include "bioplib/angle.h"
#include "bioplib/macros.h"

#include "cgutil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MODE_TARGET 0
#define MODE_ANGLE  1

#define DEF_MODE   MODE_TARGET
#define DEF_TARGET 200
#define DEF_ANGMAX (60.0 * PI / 180.0)
#define DEF_ANGLE  (30.0 * PI / 180.0)

#define ANGSTEP    (5.0 * PI / 180.0)

typedef struct _torsions
{
   struct _torsions *next;
   REAL *phi,
        *psi,
        *omega;
   int  NPhi,
        NPsi,
        NOmega;
   BOOL uniq;
}  TORSIONS;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL DoProcessing(FILE *cgin, FILE *cgout, PDB *pdb, HEADER header,
                  int mode, REAL ang, REAL MaxAngle, int target);
BOOL WriteKeptConfs(FILE *cgin, FILE *cgout, TORSIONS *torsions);
BOOL ParseCmdLine(int argc, char **argv, char *cginfile, char *pdbfile, 
                  char *cgoutfile, int *mode, REAL *ang, REAL *angmax,
                  int *target);
void Usage(void);
void PatchPDB(PDB *pdb, COOR *coor, HEADER header);
BOOL ClusterTorsions(TORSIONS *torsions, int mode, REAL ang, 
                     REAL MaxAngle, int target);
TORSIONS *CalcTorsions(PDB *pdb, HEADER header);
void DoCalcTors(PDB *start, PDB *stop, TORSIONS *torsion);
int DoCluster(TORSIONS *torsions, REAL ang);
BOOL DifferentConf(TORSIONS *t1, TORSIONS *t2, REAL ang);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program to perform clustering of CSearch conformations.

   28.11.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   char   cginfile[MAXBUFF],
          cgoutfile[MAXBUFF],
          pdbfile[MAXBUFF];
   int    mode,
          target,
          natoms,
          error=0;
   REAL   ang,
          angmax;
   FILE   *cgin,
          *cgout,
          *pdbfp;
   HEADER header;
   PDB    *pdb;
   
   if(ParseCmdLine(argc, argv, cginfile, pdbfile, cgoutfile, &mode,
                   &ang, &angmax, &target))
   {
      /* Open the files                                                 */
      if((cgin = fopen(cginfile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input CG file: %s\n",cginfile);
         return(1);
      }
      if((cgout = fopen(cgoutfile,"w"))==NULL)
      {
         fprintf(stderr,"Unable to open output CG file: %s\n",cgoutfile);
         return(1);
      }
      if((pdbfp = fopen(pdbfile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open PDB file: %s\n",pdbfile);
         return(1);
      }

      if(ReadHeader(cgin, &header))
      {
         if((pdb = ReadPDB(pdbfp, &natoms)) != NULL)
         {
            if(!DoProcessing(cgin, cgout, pdb, header, mode, ang, 
                             angmax, target))
            {
               fprintf(stderr,"Unable to complete processing\n");
               error = 1;
            }
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file\n");
            error = 1;
         }

         KillHeader(&header);
      }
      else
      {
         fprintf(stderr,"Unable to read header from CG file\n");
         error = 1;
      }
   }
   else
   {
      Usage();
   }
   
   return(error);
}


/************************************************************************/
/*>BOOL DoProcessing(FILE *cgin, FILE *cgout, PDB *pdb, HEADER header,
                     int mode, REAL ang, REAL MaxAngle, int target)
   -------------------------------------------------------------------
   Main processing section. Reads conformations, calls routines to
   do clustering and to write out the retained conformations.

   28.11.94 Original    By: ACRM
   30.11.94 Skips the reference conformation
*/
BOOL DoProcessing(FILE *cgin, FILE *cgout, PDB *pdb, HEADER header,
                  int mode, REAL ang, REAL MaxAngle, int target)
{
   int      NConf     = 0;
   TORSIONS *torsions = NULL,
            *t;
   COOR     *coor;
   REAL     energy;
   
   /* Read in all conformations and store torsion info                  */
   for(;;)
   {
      /* Read conformation and check it's OK                            */
      if((coor = GetConf(cgin, header, &energy))==NULL)
         break;
      if(coor == (COOR *)(-1))
      {
         fprintf(stderr, "Error reading conformation %d\n",NConf);
         if(NConf==0)return(FALSE);
         NConf--;
         break;
      }
      
      if(NConf++ != 0)
      {
         /* Patch the conformation into the PDB                         */
         PatchPDB(pdb, coor, header);
         
         /* Calculate the backbone torsions and store them              */
         if((torsions = CalcTorsions(pdb, header))==NULL)
         {
            fprintf(stderr, "Unable to store torsions\n");
            return(FALSE);
         }
      }

      free(coor);
   }

   /* Output number of conformations                                    */
   printf("Before clustering, conformation file contained %d confs.\n",
          NConf);

   /* Cluster the torsions                                              */
   if(!ClusterTorsions(torsions, mode, ang, MaxAngle, target))
   {
      fprintf(stderr, "Unable to cluster conformations\n");
      return(FALSE);
   }

   /* Write the confs we're keeping                                     */
   WriteKeptConfs(cgin, cgout, torsions);
}


/************************************************************************/
/*>BOOL WriteKeptConfs(FILE *cgin, FILE *cgout, TORSIONS *torsions)
   ----------------------------------------------------------------
   Write out a CG file containing unique conformations

   28.11.94 Original    By: ACRM
*/
BOOL WriteKeptConfs(FILE *cgin, FILE *cgout, TORSIONS *torsions)
{
   HEADER   header;
   TORSIONS *t;
   COOR     *coor;
   REAL     energy;
   int      NConf = 0;
   
   rewind(cgin);
   
   /* Read and write the header                                         */
   if(ReadHeader(cgin, &header))
   {
      WriteHeader(cgout, header);

      /* Read and write the reference conformation                      */
      if((coor = GetConf(cgin, header, &energy))!=NULL)
      {
         if(coor == (COOR *)(-1))
         {
            fprintf(stderr, "Error reading conformation %d\n",NConf);
            return(FALSE);
         }
         NConf++;
         
         PutConf(cgout, header, coor, energy);
         free(coor);

         /* Step through the torsion linked list, reading the equivalent
            conformations and writing them if the unique flag is set
         */
         for(t=torsions;;NEXT(t),NConf++)
         {
            /* Read conformation and check it's OK                      */
            if((coor = GetConf(cgin, header, &energy))==NULL)
               break;
            if(coor == (COOR *)(-1))
            {
               /* We don't issue an error message, since we've already
                  done so on the first pass
               */
               break;
            }
            
            if(t->uniq)
               PutConf(cgout, header, coor, energy);
            
            free(coor);
         }
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *cginfile, 
                     char *pdbfile, char *cgoutfile, int *mode, 
                     REAL *ang, REAL *angmax, int *target)
   ------------------------------------------------------------
   Parse the command line

   28.11.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *cginfile, char *pdbfile, 
                  char *cgoutfile, int *mode, REAL *ang, REAL *angmax,
                  int *target)
{
   argc--;
   argv++;
   
   *mode        = DEF_MODE;
   *target      = DEF_TARGET;
   *angmax      = DEF_ANGMAX;
   *ang         = DEF_ANGLE;

   cginfile[0]  = '\0';
   cgoutfile[0] = '\0';
   pdbfile[0]   = '\0';

   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
            argc--;
            argv++;
            if(!sscanf(argv[0],"%d",target))
               return(FALSE);
            *mode = MODE_TARGET;
            break;
         case 'a':
            argc--;
            argv++;
            if(!sscanf(argv[0],"%lf",ang))
               return(FALSE);
            *ang  *= PI/180.0;
            *mode  = MODE_ANGLE;
            break;
         case 'm':
            argc--;
            argv++;
            if(!sscanf(argv[0],"%lf",angmax))
               return(FALSE);
            *angmax *= PI/180.0;
            *mode   = MODE_TARGET;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 3 arguments left                       */
         if(argc != 3)
            return(FALSE);
         
         /* Copy the filenames                                          */
         strcpy(cginfile,  argv[0]);
         strcpy(pdbfile,   argv[1]);
         strcpy(cgoutfile, argv[2]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message.

   28.11.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nClusterCG V1.1 (c) Dr. Andrew C.R. Martin, UCL\n\n");
   
   fprintf(stderr,"clustercg [-a angle] [-n target] [-m maxang] in.cg \
ref.pdb out.cg\n");
   fprintf(stderr,"          -a Specify initial angular cutoff to use  \
(Default: %.1f)\n",180.0*DEF_ANGLE/PI);
   fprintf(stderr,"          -n Specify target number of \
conformations (Default: %d)\n",DEF_TARGET);
   fprintf(stderr,"          -m Specify max angular cutoff for -n      \
(Default: %.1f)\n\n",180.0*DEF_ANGMAX/PI);

   fprintf(stderr,"ClusterCG rejects conformations from a CSearch CG \
file which are similar\n");
   fprintf(stderr,"to other conformations on the basis of similarity in \
phi, psi and\n");
   fprintf(stderr,"omega torsion angles. By default, it will gradually \
increase the cutoff\n");
   fprintf(stderr,"until a target number of conformations is left.\n");
   fprintf(stderr,"   The -a option will cause a single round of \
clustering to be performed\n");
   fprintf(stderr,"using the specified angle. If -a is FOLLOWED by -n, \
then the angle\n");
   fprintf(stderr,"specified for -a will be the initial angle used for \
clustering. The\n");
   fprintf(stderr,"-m option specifies the maximum allowed angle while \
trying to reach a\n");
   fprintf(stderr,"target number of conformations.\n\n");
}


/************************************************************************/
/*>void PatchPDB(PDB *pdb, COOR *coor, HEADER header)
   --------------------------------------------------
   Runs through a PDB linked list patching in atoms from the coor array
   28.11.94 Original    By: ACRM
*/
void PatchPDB(PDB *pdb, COOR *coor, HEADER header)
{
   int i, atomnum;
   PDB *p;

   for(p=pdb, atomnum=1; p!=NULL; NEXT(p), atomnum++)
   {
      for(i=0; i<header.ncons; i++)
      {
         if(header.ConsAtoms[i] == atomnum)
         {
            p->x = coor[i].x;
            p->y = coor[i].y;
            p->z = coor[i].z;
            break;
         }
      }
   }
}


/************************************************************************/
/*>TORSIONS *CalcTorsions(PDB *pdb, HEADER header)
   -----------------------------------------------
   Run through a conformation, create storage space and call the routine
   to calculate and store the torsions.

   28.11.94 Original    By: ACRM
*/
TORSIONS *CalcTorsions(PDB *pdb, HEADER header)
{
   static TORSIONS *torsions = NULL,
                   *t;
   static int      NCAlpha = 0,
                   MinAtom,
                   MaxAtom;
   int             natom,
                   i;
   PDB             *start, 
                   *stop;
   static PDB      **indx = NULL;
   
   /* Allocate space in the linked list                                 */
   if(torsions==NULL)
   {
      INIT(torsions, TORSIONS);
      t = torsions;

      /* This is the first call, so find the number of CAs in the
         constructed region
      */
      if((indx = IndexPDB(pdb, &natom))==NULL)
         return(NULL);

      MinAtom = MaxAtom = header.ConsAtoms[0];
      
      for(i=0; i<header.ncons; i++)
      {
         if(!strncmp((indx[header.ConsAtoms[i]])->atnam, "CA  ", 4))
         {
            NCAlpha++;
         }
         if(header.ConsAtoms[i] < MinAtom) 
            MinAtom = header.ConsAtoms[i];
         if(header.ConsAtoms[i] > MaxAtom) 
            MaxAtom = header.ConsAtoms[i];
      }
   }
   else
   {
      ALLOCNEXT(t, TORSIONS);
   }
   if(t==NULL) return(NULL);
   
   /* Allocate space for torsions                                       */
   if((t->phi = (REAL *)malloc(NCAlpha * sizeof(REAL)))==NULL)
      return(NULL);
   if((t->psi = (REAL *)malloc(NCAlpha * sizeof(REAL)))==NULL)
      return(NULL);
   if((t->omega = (REAL *)malloc(NCAlpha * sizeof(REAL)))==NULL)
      return(NULL);

   t->NPhi = t->NPsi = t->NOmega = 0;
   t->uniq = TRUE;

   /* Run between min and max                                           */
   start = indx[MinAtom];
   stop  = indx[MaxAtom];

   DoCalcTors(start, stop, t);

   return(torsions);
}


/************************************************************************/
/*>void DoCalcTors(PDB *start, PDB *stop, TORSIONS *torsion)
   ---------------------------------------------------------
   Calculate and store the torsion angles for a conformation

   28.11.94 Original    By: ACRM
*/
void DoCalcTors(PDB *start, PDB *stop, TORSIONS *torsion)
{
   PDB *p,
       *Nprev  = NULL,
       *CAprev = NULL,
       *Cprev  = NULL,
       *N      = NULL,
       *CA     = NULL,
       *C      = NULL;
   
   for(p=start; p!=stop->next; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4))
      {
         Nprev = N;
         N = p;
         if(Nprev != NULL && CAprev != NULL && Cprev != NULL)
         {
            torsion->psi[(torsion->NPsi)++] = 
               phi(Nprev->x,  Nprev->y,  Nprev->z,
                   CAprev->x, CAprev->y, CAprev->z,
                   Cprev->x,  Cprev->y,  Cprev->z,
                   N->x,      N->y,      N->z);
         }
      }
      else if(!strncmp(p->atnam,"CA  ",4))
      {
         CAprev = CA;
         CA = p;
         if(CAprev != NULL && Cprev != NULL && N != NULL)
         {
            torsion->omega[(torsion->NOmega)++] = 
               phi(CAprev->x, CAprev->y, CAprev->z,
                   Cprev->x,  Cprev->y,  Cprev->z,
                   N->x,      N->y,      N->z,
                   CA->x,     CA->y,     CA->z);
         }
      }
      else if(!strncmp(p->atnam,"C   ",4))
      {
         Cprev = C;
         C = p;
         if(Cprev != NULL && N != NULL && CA != NULL)
         {
            torsion->phi[(torsion->NPhi)++] = 
               phi(Cprev->x,  Cprev->y,  Cprev->z,
                   N->x,      N->y,      N->z,
                   CA->x,     CA->y,     CA->z,
                   C->x,      C->y,      C->z);
         }
      }
   }
}


/************************************************************************/
/*>BOOL ClusterTorsions(TORSIONS *torsions, int mode, REAL InitAngle, 
                        REAL MaxAngle, int target)
   ------------------------------------------------------------------
   Performs the main function of rejecting similar conformations.
   If mode is MODE_ANGLE, performs a single pass using InitAngle as the
   resolution cutoff. If mode is MODE_TARGET, performs repeated
   iterations increasing the cutoff from angle to MaxAngle until the
   target number of conformations is reached.

   30.11.94 Original    By: ACRM
*/
BOOL ClusterTorsions(TORSIONS *torsions, int mode, REAL InitAngle, 
                     REAL MaxAngle, int target)
{
   int  NKept;
   REAL ang;
   
   switch(mode)
   {
   case MODE_ANGLE:
      NKept = DoCluster(torsions, InitAngle);
      printf("Cluster at %.1f degrees. Retained %d conformations.\n",
             180.0 * InitAngle / PI, NKept);
      break;
   case MODE_TARGET:
      for(ang=InitAngle; ang<=MaxAngle; ang += ANGSTEP)
      {
         NKept = DoCluster(torsions, ang);
         printf("Cluster at %f degrees. Retained %d conformations.\n",
                180.0 * ang / PI, NKept);
         if(NKept <= target)
            break;
      }

      if(NKept > target)
         printf("Target not reached, maximum angular cutoff exceeded.\n");

      break;
   default:
      break;
   }
      
   return(TRUE);
}


/************************************************************************/
/*>int DoCluster(TORSIONS *torsions, REAL ang)
   -------------------------------------------
   Does the actual clustering work. This is not a true clustering since
   we simply look for similar conformations and select the first one
   of a batch which we find.

   30.11.94 Original   By: ACRM
*/
int DoCluster(TORSIONS *torsions, REAL ang)
{
   TORSIONS *t1,
            *t2;
   int      NUnique;
   
   /* Walk through the conformations                                    */
   for(t1=torsions; t1!=NULL; NEXT(t1))
   {
      /* If this conf is unique                                         */
      if(t1->uniq)
      {
         /* Walk through the rest of the conformations                  */
         for(t2=t1->next; t2!=NULL; NEXT(t2))
         {
            /* If conf t2 is currently unique, see if it is similar to
               the current t1
            */
            if(t2->uniq)
               t2->uniq = DifferentConf(t1, t2, ang);
         }
      }
   }

   /* Count the number of unique conformations remaining                */
   for(t1=torsions, NUnique = 0; t1!=NULL; NEXT(t1))
   {
      if(t1->uniq)
         NUnique++;
   }

   return(NUnique);
}


/************************************************************************/
/*>BOOL DifferentConf(TORSIONS *t1, TORSIONS *t2, REAL ang)
   --------------------------------------------------------
   Determines whether two conformations are the same within the specified
   angular cutoff

   30.11.94 Original    By: ACRM
   01.02.95 Corrected pi/2 and pi to pi and 2.0 * pi ...
*/
BOOL DifferentConf(TORSIONS *t1, TORSIONS *t2, REAL ang)
{
   int  i;
   REAL TorDiff,
        pi = PI;
   
   /* Do the Phi angles                                                 */
   for(i=0; i<(t1->NPhi); i++)
   {
      TorDiff = ABS((t1->phi[i]) - (t2->phi[i]));

      if(TorDiff > pi) 
         TorDiff = (2.0 * pi) - TorDiff;

      if(TorDiff > ang)
         return(TRUE);
   }

   /* Do the Psi angles                                                 */
   for(i=0; i<(t1->NPsi); i++)
   {
      TorDiff = ABS((t1->psi[i]) - (t2->psi[i]));

      if(TorDiff > pi) 
         TorDiff = (2.0 * pi) - TorDiff;

      if(TorDiff > ang)
         return(TRUE);
   }

   /* Do the Omega angles                                               */
   for(i=0; i<(t1->NOmega); i++)
   {
      TorDiff = ABS((t1->omega[i]) - (t2->omega[i]));

      if(TorDiff > pi) 
         TorDiff = (2.0 * pi) - TorDiff;

      if(TorDiff > ang)
         return(TRUE);
   }

   /* All torsions were within prescribed cutoff, so conformations are not
      different
   */
   return(FALSE);
}
