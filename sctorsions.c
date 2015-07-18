/*************************************************************************

   Program:    sctorsions
   File:       sctorsions.c
   
   Version:    V1.0
   Date:       10.10.01
   Function:   Take a CF-CONGEN format .CG file and a coordinate file
               and print torsions and energies
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2001
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
   EMail:      andrew@bioinf.org.uk
               
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

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/fsscanf.h"
#include "bioplib/angle.h"
#include "cgutil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256
#define MAXRESTYPE 20

typedef struct _torsions
{
   struct _torsions *next;
   PDB *atm1,
       *atm2,
       *atm3,
       *atm4;
} TORSIONS;

/************************************************************************/
/* Globals
*/
PDB **gPointers = NULL;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL RunSCTorsions(FILE *cg_fp, FILE *crd_fp, FILE *out);
BOOL TestEachConformation(FILE *out, FILE *cg_fp, PDB *pdb, 
                          HEADER header);
BOOL PatchAConformation(PDB *pdb, HEADER header, COOR *coor);
BOOL PrintTorsionsAndEnergies(FILE *out, PDB *pdb, HEADER header, 
                              REAL energy);
PDB *ReadCrd(FILE *fp, int *ncoor);
BOOL ParseCmdLine(int argc, char **argv, char *cgfile, char *crdfile, 
                  char *outfile);
void Usage(void);
BOOL CreateTorAtomList(char *resnam, char toratomlist[5][4][8]);
char *PrintResidue(PDB *pdb);
char *PrintAtomQuartet(char *atm1, char *atm2, char *atm3, char *atm4);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Input:     
   Output:    
   Returns:   

   10.10.01 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *cg_fp  = NULL,
        *crd_fp = NULL,
        *out    = stdout;
   char cgfile[MAXBUFF],
        crdfile[MAXBUFF],
        outfile[MAXBUFF];

   cgfile[0] = crdfile[0] = outfile[0] = '\0';
   
   if(ParseCmdLine(argc, argv, cgfile, crdfile, outfile))
   {
      if(cgfile[0] && crdfile[0])
      {
         if((cg_fp = fopen(cgfile,"r"))!=NULL)
         {
            if((crd_fp = fopen(crdfile,"r"))!=NULL)
            {
               if(outfile[0])
               {
                  if((out=fopen(outfile,"w"))==NULL)
                  {
                     fprintf(stderr,"Can't write %s\n", outfile);
                     return(1);
                  }
               }
               /* All files opened OK, so run the main program          */

               return(RunSCTorsions(cg_fp, crd_fp, out));
            }
            else
            {
               fprintf(stderr,"Can't open %s\n", crdfile);
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"Can't open %s\n", cgfile);
            return(1);
         }
      }
      else
      {
         Usage();
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL RunSCTorsions(FILE *cg_fp, FILE *crd_fp, FILE *out)
   --------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Main routine that extracts reads the header and coordinate data and
   then calls TestEachConformation() to process each conformation in turn

   10.10.01 Original   By: ACRM
*/
BOOL RunSCTorsions(FILE *cg_fp, FILE *crd_fp, FILE *out)
{
   HEADER header;
   BOOL   retval = FALSE;
   PDB    *pdb = NULL;
   int    ncoor;
   
   if(ReadHeader(cg_fp, &header))
   {
      if((pdb = ReadCrd(crd_fp, &ncoor))!=NULL)
      {
         retval = TestEachConformation(out, cg_fp, pdb, header);
         FREELIST(pdb, PDB);
      }
      else
      {
         fprintf(stderr,"Error reading CHARMM .CRD file\n");
      }
   }
   else
   {
      fprintf(stderr,"Error reading CG header\n");
   }
   
   return(retval);
}

   
/************************************************************************/
/*>BOOL TestEachConformation(FILE *out, FILE *cg_fp, PDB *pdb, 
                             HEADER header)
   -----------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Reads each conformation from a CG file and patches the coordinates
   then calls PrintTorsionsAndEnergies() to calculate torsions and
   energy

   10.10.01 Original   By: ACRM
*/
BOOL TestEachConformation(FILE *out, FILE *cg_fp, PDB *pdb, HEADER header)
{
   COOR *coor;
   REAL energy;
   int  count = 0;
      
   while((coor=GetConf(cg_fp, header, &energy))!=NULL)
   {
      if(coor == (COOR *)(-1))
      {
         fprintf(stderr,"Error in reading coordinates from CG file\n");
         return(FALSE);
      }
      
      /* We have a new malloc'd coordinate set in coor                  */
      if(!PatchAConformation(pdb, header, coor))
      {
         return(FALSE);
      }

      fprintf(out, "Conf %d : ", count++);
      if(!PrintTorsionsAndEnergies(out, pdb, header, energy))
      {
         return(FALSE);
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL PatchAConformation(PDB *pdb, HEADER header, COOR *coor)
   ------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Patches the coordinate file with coordinates from a given conformation

   10.10.01 Original   By: ACRM
*/
BOOL PatchAConformation(PDB *pdb, HEADER header, COOR *coor)
{
   PDB *p;
   int i;
   
   /* If it's the first call to this routine we find pointers into the PDB
      linked list for the atoms
   */
   if(gPointers == NULL)
   {
      BOOL *found;
      if((found = (BOOL *)malloc(header.ncons * sizeof(BOOL)))==NULL)
      {
         fprintf(stderr,"No memory for finding constructed atoms\n");
         return(FALSE);
      }
      for(i=0; i<header.ncons; i++)
      {
         found[i] = FALSE;
      }
      
      if((gPointers = (PDB **)malloc(header.ncons * sizeof(PDB *)))==NULL)
      {
         fprintf(stderr,"No memory for constructed atom pointers\n");
         return(FALSE);
      }
      for(p=pdb; p!=NULL; NEXT(p))
      {
         for(i=0; i<header.ncons; i++)
         {
            if(p->atnum == header.ConsAtoms[i])
            {
               gPointers[i] = p;
               found[i] = TRUE;
            }
         }
      }

      /* Check we found all the atoms                                   */
      for(i=0; i<header.ncons; i++)
      {
         if(!found[i])
         {
            fprintf(stderr,"Atom %d not found in coordinate file\n",
                    header.ConsAtoms[i]);
            return(FALSE);
         }
      }
      free(found);
   }

   for(i=0; i<header.ncons; i++)
   {
      gPointers[i]->x = coor[i].x;
      gPointers[i]->y = coor[i].y;
      gPointers[i]->z = coor[i].z;
   }

   return(TRUE);
}


/************************************************************************/
/*>PDB *ReadCrd(FILE *fp, int *ncoor)
   ----------------------------------
   Input:     
   Output:    
   Returns:   

   Reads a coordinate set from a CRD file into a PDB linked list

   10.10.01 Original   By: ACRM
*/
PDB *ReadCrd(FILE *fp, int *ncoor)
{
   char buffer[MAXBUFF];
   PDB  *pdb = NULL, 
        *p;
   
   /* Read and discard header                                           */
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      KILLTRAILSPACES(buffer);
      if(!strcmp(buffer,"*"))
         break;
   }
   
   /* Read the number of coordinates                                    */
   if(!fgets(buffer, MAXBUFF, fp))
      return(NULL);
   TERMINATE(buffer);
   if(!sscanf(buffer, "%d", ncoor))
      return(NULL);
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      if(pdb==NULL)
      {
         INIT(pdb, PDB);
         p = pdb;
      }
      else
      {
         ALLOCNEXT(p, PDB);
      }

      if(p==NULL)
      {
         fprintf(stderr,"No memory for PDB linked list\n");
         FREELIST(pdb, PDB);
         return(NULL);
      }
      
      fsscanf(buffer, "%5d%5d%1x%4s%1x%4s%10lf%10lf%10lf",
              &(p->atnum), &(p->resnum), p->resnam,
              p->atnam, &(p->x), &(p->y), &(p->z));
      p->occ = p->bval = (REAL)0.0;
      strcpy(p->junk, "ATOM  ");
      strcpy(p->atnam_raw, p->atnam);
      strcpy(p->insert, " ");
      strcpy(p->chain, " ");
   }
   
   return(pdb);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Input:     
   Output:    
   Returns:   

   Prints a Usage message

   10.10.01 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nsctorsions V1.0 (c) Dr. Andrew C.R. Martin, \
University of Reading\n");

   fprintf(stderr,"\nUsage: sctorsions cgfile crdfile [output]\n");
   fprintf(stderr,"       cgfile  - a conformation file in CHARMM-free \
congen ASCII format\n");
   fprintf(stderr,"       crdfile - a CHARMM CRD format coordinate file \
with all hydrogens\n");
   fprintf(stderr,"                 NTER and CTER residues added (from \
a 'genpsf' run of \n");
   fprintf(stderr,"                 CONGEN)\n");
   fprintf(stderr,"       output  - an optional output file (output is \
to stdout if not\n");
   fprintf(stderr,"                 specified)\n");

   fprintf(stderr,"\nsctorsions reads a CHARMM-free congen ASCII format \
CG file (generated \n");
   fprintf(stderr,"from a CONGEN CG file using cg2cs) and produces an \
output file listing\n");
   fprintf(stderr,"the residues involved, their torsion angles and the \
energy calculated\n");
   fprintf(stderr,"by CONGEN. Note that many energies may be generated \
for each combination\n");
   fprintf(stderr,"of similar chi torsion angles. This is because the \
hydrogen positions\n");
   fprintf(stderr,"may also be rotated (e.g. in a serine) and the \
angles for these \n");
   fprintf(stderr,"H-torsions are not displayed. Simple rule is just \
take the lowest\n");
   fprintf(stderr,"energy from each group!\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *cgfile, char *crdfile, 
                     char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *cgfile      Name of CG file
            char   *crdfile     Name of .CRD file
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   10.10.01 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *cgfile, char *crdfile, 
                  char *outfile)
{
   argc--;
   argv++;

   cgfile[0] = crdfile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 2 or 3 arguments left             */
         if((argc > 3) || (argc < 2))
            return(FALSE);
         
         /* Copy the first to cgfile and the second to crdfile          */
         strcpy(cgfile, argv[0]);
         argc--;
         argv++;
         strcpy(crdfile, argv[0]);
         argc--;
         argv++;
         
         /* If there's another, copy it to outfile                      */
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL PrintTorsionsAndEnergies(FILE *out, PDB *pdb, HEADER header, 
                                 REAL energy)
   ------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Determine the torsion angles for the residues containing each of the
   atoms. Print them and the associated energy.

   10.10.01 Original   By: ACRM
*/
BOOL PrintTorsionsAndEnergies(FILE *out, PDB *pdb, HEADER header, 
                              REAL energy)
{
   PDB *p,
       *firstAtom,
       *atms[8];
   BOOL *done;
   int i,j;
   REAL angle;
   TORSIONS *torsions = NULL;
   char toratomlist[5][4][8];

   if(gPointers == NULL)
   {
      fprintf(stderr,"List of reconstructed atoms pointers has not been \
built\n");
      return(FALSE);
   }
   
   /* Create an array to flag whether each atom has been assigned       */
   if((done = (BOOL *)malloc(header.ncons * sizeof(BOOL)))==NULL)
   {
      fprintf(stderr,"No memory for atom/torsion assignment flags\n");
      return(FALSE);
   }
   for(i=0;i<header.ncons; i++)
   {
      done[i] = FALSE;
   }
   
   /* Run through the spun atoms                                        */
   for(i=0; i<header.ncons; i++)
   {
      /* If this atom not yet assigned to a torsion                     */
      if(!done[i])
      {
         /* Run through the PDB list to find the same residue           */
         firstAtom = NULL;
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if((p->resnum == gPointers[i]->resnum) &&
               (p->chain[0] == gPointers[i]->chain[0]) &&
               (p->insert[0] == gPointers[i]->chain[0]))
            {
               firstAtom = p;
               break;
            }
         }
         if(firstAtom==NULL)
         {
            fprintf(stderr,"Internal inconsistency! Can't find residue \
for constructed atoms\n");
            FREELIST(torsions, TORSIONS);
            return(FALSE);
         }
         
         /* Set the done flags to mark all atoms in this residue as being
            done for torsions 
         */
         for(j=0; j<header.ncons; j++)
         {
            if(!done[j])
            {
               if((gPointers[j]->resnum == gPointers[i]->resnum) &&
                  (gPointers[j]->chain[0] == gPointers[j]->chain[0]) &&
                  (gPointers[j]->insert[0] == gPointers[j]->insert[0]))
               {
                  done[j] = TRUE;
               }
            }
         }

         /* Find the list of atoms involved in all torsions             */
         CreateTorAtomList(firstAtom->resnam, toratomlist);

         /* For each possible torsion angle                             */
         for(i=0; i<5; i++)
         {
            /* Break out if there are no more torsions                  */
            if(toratomlist[i][0][0] == ' ')
               break;
            
            /* Find pointers to the atoms                               */
            for(j=0; j<4; j++)
            {
               atms[j] = FindAtomInRes(firstAtom, toratomlist[i][j]);
               if(atms[j] == NULL)
               {
                  fprintf(stderr,"Expected torsion atom, %s, not found \
in %s\n",
                          toratomlist[i][j], PrintResidue(firstAtom));
               }
            }
            
            /* Calculate the torsion angle                              */
            angle = phi(atms[0]->x, atms[0]->y, atms[0]->z,
                        atms[1]->x, atms[1]->y, atms[1]->z,
                        atms[2]->x, atms[2]->y, atms[2]->z,
                        atms[3]->x, atms[3]->y, atms[3]->z);
            
            /* Print the result                                         */
            fprintf(stdout, "%s Chi %d (%s) = %.4f : ",
                    PrintResidue(firstAtom),
                    i+1, 
                    PrintAtomQuartet(toratomlist[i][0], 
                                     toratomlist[i][1], 
                                     toratomlist[i][2], 
                                     toratomlist[i][3]),
                    (REAL)180.0 * angle / PI);
            
         }  /* End of run through all possible chis                     */
      }  /* This atom's residue has not yet been done                   */
   }  /* End of run through all constructed atoms                       */
   
   fprintf(stdout, "%f\n", energy);
      
   /* Free temporary storage                                            */
   free(done);
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL CreateTorAtomList(char *resnam, char toratomlist[5][4][8])
   ---------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Contains lookup tables for the atoms involved in torsions in each of
   the amino acid types. Outputs an array containing the atom lists
   for the required type.

   10.10.01 Original   By: ACRM
*/
BOOL CreateTorAtomList(char *resnam, char toratomlist[5][4][8])
{
   int i, j;
   
   
   /* Initialise lookup tables                                          */
   static char residues[MAXRESTYPE][8] = 
   {"ALA ","CYS ","ASP ","GLU ","PHE ",
    "GLY ","HIS ","ILE ","LYS ","LEU ",
    "MET ","ASN ","PRO ","GLN ","ARG ",
    "SER ","THR ","VAL ","TRP ","TYR "
   };
   static char chi1[MAXRESTYPE][4][8] = 
   {{"    ", "    ", "    ", "    "}, /* ALA */
    {"N   ", "CA  ", "CB  ", "SG  "}, /* CYS */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* ASP */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* GLU */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* PHE */
    {"    ", "    ", "    ", "    "}, /* GLY */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* HIS */
    {"N   ", "CA  ", "CB  ", "CG1 "}, /* ILE */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* LYS */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* LEU */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* MET */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* ASN */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* PRO */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* GLN */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* ARG */
    {"N   ", "CA  ", "CB  ", "OG  "}, /* SER */
    {"N   ", "CA  ", "CB  ", "OG1 "}, /* THR */
    {"N   ", "CA  ", "CB  ", "CG1 "}, /* VAL */
    {"N   ", "CA  ", "CB  ", "CG  "}, /* TRP */
    {"N   ", "CA  ", "CB  ", "CG  "}  /* TYR */
   };
   static char chi2[MAXRESTYPE][4][8] = 
   {{"    ", "    ", "    ", "    "}, /* ALA */
    {"    ", "    ", "    ", "    "}, /* CYS */
    {"CA  ", "CB  ", "CG  ", "OD1 "}, /* ASP */
    {"CA  ", "CB  ", "CG  ", "CD  "}, /* GLU */
    {"CA  ", "CB  ", "CG  ", "CD1 "}, /* PHE */
    {"    ", "    ", "    ", "    "}, /* GLY */
    {"CA  ", "CB  ", "CG  ", "ND1 "}, /* HIS */
    {"CA  ", "CB  ", "CG1 ", "CD  "}, /* ILE Note CD not CD1 (CHARMM!)  */
    {"CA  ", "CB  ", "CG  ", "CD  "}, /* LYS */
    {"CA  ", "CB  ", "CG  ", "CD1 "}, /* LEU */
    {"CA  ", "CB  ", "CG  ", "SD  "}, /* MET */
    {"CA  ", "CB  ", "CG  ", "OD1 "}, /* ASN */
    {"CA  ", "CB  ", "CG  ", "CD  "}, /* PRO */
    {"CA  ", "CB  ", "CG  ", "CD  "}, /* GLN */
    {"CA  ", "CB  ", "CG  ", "CD  "}, /* ARG */
    {"    ", "    ", "    ", "    "}, /* SER */
    {"    ", "    ", "    ", "    "}, /* THR */
    {"    ", "    ", "    ", "    "}, /* VAL */
    {"CA  ", "CB  ", "CG  ", "CD1 "}, /* TRP */
    {"CA  ", "CB  ", "CG  ", "CD1 "}  /* TYR */
   };
   static char chi3[MAXRESTYPE][4][8] = 
   {{"    ", "    ", "    ", "    "}, /* ALA */
    {"    ", "    ", "    ", "    "}, /* CYS */
    {"    ", "    ", "    ", "    "}, /* ASP */
    {"CB  ", "CG  ", "CD  ", "OE1 "}, /* GLU */
    {"    ", "    ", "    ", "    "}, /* PHE */
    {"    ", "    ", "    ", "    "}, /* GLY */
    {"    ", "    ", "    ", "    "}, /* HIS */
    {"    ", "    ", "    ", "    "}, /* ILE */
    {"CB  ", "CG  ", "CD  ", "CE  "}, /* LYS */
    {"    ", "    ", "    ", "    "}, /* LEU */
    {"CB  ", "CG  ", "SD  ", "CE  "}, /* MET */
    {"    ", "    ", "    ", "    "}, /* ASN */
    {"    ", "    ", "    ", "    "}, /* PRO */
    {"CB  ", "CG  ", "CD  ", "OE1 "}, /* GLN */
    {"CB  ", "CG  ", "CD  ", "NE  "}, /* ARG */
    {"    ", "    ", "    ", "    "}, /* SER */
    {"    ", "    ", "    ", "    "}, /* THR */
    {"    ", "    ", "    ", "    "}, /* VAL */
    {"    ", "    ", "    ", "    "}, /* TRP */
    {"    ", "    ", "    ", "    "}  /* TYR */
   };
   static char chi4[MAXRESTYPE][4][8] = 
   {{"    ", "    ", "    ", "    "}, /* ALA */
    {"    ", "    ", "    ", "    "}, /* CYS */
    {"    ", "    ", "    ", "    "}, /* ASP */
    {"    ", "    ", "    ", "    "}, /* GLU */
    {"    ", "    ", "    ", "    "}, /* PHE */
    {"    ", "    ", "    ", "    "}, /* GLY */
    {"    ", "    ", "    ", "    "}, /* HIS */
    {"    ", "    ", "    ", "    "}, /* ILE */
    {"CG  ", "CD  ", "CE  ", "NZ  "}, /* LYS */
    {"    ", "    ", "    ", "    "}, /* LEU */
    {"    ", "    ", "    ", "    "}, /* MET */
    {"    ", "    ", "    ", "    "}, /* ASN */
    {"    ", "    ", "    ", "    "}, /* PRO */
    {"    ", "    ", "    ", "    "}, /* GLN */
    {"CG  ", "CD  ", "NE  ", "CZ  "}, /* ARG */
    {"    ", "    ", "    ", "    "}, /* SER */
    {"    ", "    ", "    ", "    "}, /* THR */
    {"    ", "    ", "    ", "    "}, /* VAL */
    {"    ", "    ", "    ", "    "}, /* TRP */
    {"    ", "    ", "    ", "    "}  /* TYR */
   };
   static char chi5[MAXRESTYPE][4][8] = 
   {{"    ", "    ", "    ", "    "}, /* ALA */
    {"    ", "    ", "    ", "    "}, /* CYS */
    {"    ", "    ", "    ", "    "}, /* ASP */
    {"    ", "    ", "    ", "    "}, /* GLU */
    {"    ", "    ", "    ", "    "}, /* PHE */
    {"    ", "    ", "    ", "    "}, /* GLY */
    {"    ", "    ", "    ", "    "}, /* HIS */
    {"    ", "    ", "    ", "    "}, /* ILE */
    {"    ", "    ", "    ", "    "}, /* LYS */
    {"    ", "    ", "    ", "    "}, /* LEU */
    {"    ", "    ", "    ", "    "}, /* MET */
    {"    ", "    ", "    ", "    "}, /* ASN */
    {"    ", "    ", "    ", "    "}, /* PRO */
    {"    ", "    ", "    ", "    "}, /* GLN */
    {"CD  ", "NE  ", "CZ  ", "NH1 "}, /* ARG */
    {"    ", "    ", "    ", "    "}, /* SER */
    {"    ", "    ", "    ", "    "}, /* THR */
    {"    ", "    ", "    ", "    "}, /* VAL */
    {"    ", "    ", "    ", "    "}, /* TRP */
    {"    ", "    ", "    ", "    "}  /* TYR */
   };

   for(i=0; i<MAXRESTYPE; i++)
   {
      if(!strncmp(resnam, residues[i], 3))
      {
         for(j=0; j<4; j++)
         {
            strcpy(toratomlist[0][j], chi1[i][j]);
            strcpy(toratomlist[1][j], chi2[i][j]);
            strcpy(toratomlist[2][j], chi3[i][j]);
            strcpy(toratomlist[3][j], chi4[i][j]);
            strcpy(toratomlist[4][j], chi5[i][j]);
         }
         return(TRUE);
      }
   }

   fprintf(stderr,"Torsion atom list not found for %s\n", resnam);
   return(FALSE);
}


/************************************************************************/
/*>char *PrintResidue(PDB *pdb)
   ----------------------------
   Input:     
   Output:    
   Returns:   

   Creates a string to identify a residue

   10.10.01 Original   By: ACRM
*/
char *PrintResidue(PDB *pdb)
{
   static char buffer[MAXBUFF];
   char buff2[MAXBUFF];
   
   strcpy(buffer,pdb->resnam);
   KILLTRAILSPACES(buffer);
   strcat(buffer,"-");
   strcat(buffer,pdb->chain);
   KILLTRAILSPACES(buffer);
   sprintf(buff2, "%d", pdb->resnum);
   strcat(buffer,buff2);
   KILLTRAILSPACES(buffer);
   strcat(buffer,pdb->insert);
   KILLTRAILSPACES(buffer);
   
   return(buffer);
}


/************************************************************************/
/*>char *PrintAtomQuartet(char *atm1, char *atm2, char *atm3, char *atm4)
   ----------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Creates a string concatenating 4 atom names with - signs in between

   10.10.01 Original   By: ACRM
*/
char *PrintAtomQuartet(char *atm1, char *atm2, char *atm3, char *atm4)
{
   static char buffer[MAXBUFF];
   
   strcpy(buffer, atm1);
   KILLTRAILSPACES(buffer);
   strcat(buffer, "-");

   strcat(buffer, atm2);
   KILLTRAILSPACES(buffer);
   strcat(buffer, "-");

   strcat(buffer, atm3);
   KILLTRAILSPACES(buffer);
   strcat(buffer, "-");

   strcat(buffer, atm4);
   KILLTRAILSPACES(buffer);

   return(buffer);
}
