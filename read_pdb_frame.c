#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"
#include "header.h"

int read_line(FILE *fp, int *p_ichar);

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position, int to_space );

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

double atomic_bscat_list( char *element );

int read_pdb_frame(FILE *file_fp, atom *p_molecule, int num_atoms, int *p_nstep, 
                   double *p_tstep, int pbc, double *p_latt_vec, int have_config,
                   int levcfg ) 
{

int good_read, num_of_chars, itsanum;
int ichar[LINESIZ], ndigi, sign, to_space;

int is_molecules, is_num_this_mol, is_demarcation;
int last_start, place, num, occurances, iloop, jloop;
int ic;

int keytrj, imcon, natms, dummy;

double velocity[3], force[3];

char *p_key, cchar[LINESIZ];
int ijj;

/*********************************************************************/
/**** Read in a frame from the HISTORY file **************************/
/*********************************************************************/

DEBUG=TRUE;
  if (DEBUG)
    {
      printf("Arrived in read_pdb_frame with:\n");
      printf("pbc = %d, levcfg = %d, num_atoms = %d\n", pbc, levcfg, num_atoms);
    }

 if (!have_config)
   {
     num_of_chars= read_line(file_fp, &ichar[0]);
      if (DEBUG)
        {
          printf("READ: ");
          for (iloop=0; iloop < num_of_chars; iloop++) putchar(ichar[iloop]); 
          printf("\n");
        }

     good_read= num_of_chars != -10 && locate_string( "timestep", &ichar[0], num_of_chars );   
/*     printf("GOOD_READ=%d\n", good_read); */

/*** Return if this is the last frame ********************************/

     if (!good_read) 
       {
         printf("Failed in read_hist_frame. Last line read:\n");
         for (iloop=0; iloop < num_of_chars; iloop++) putchar(ichar[iloop]); 
         printf("\nnum_of_chars = %d\n", num_of_chars);
         DEBUG=FALSE;
         return -10;
       }

/**** Current timestep number *****/
     place = 0;
     *p_nstep  = get_int(&ichar[0], &place, &itsanum,
                         &ndigi, num_of_chars, &sign);

/**** number of atoms in this frame ****/
     natms     = get_int(&ichar[0], &place, &itsanum,
                         &ndigi, num_of_chars, &sign);

     if (natms != num_atoms+1)
       {
         printf("Warning: Frame found without the correct number of atoms present\n");
         printf("         Read %d from HISTORY timestep line and expected %d from FIELD file\n",
                                    natms, num_atoms+1);

/*** Addition July 25th, If the HISTORY does not agree with FIELD believe FIELD ****/
/***                     This allows large files containing one molecule we are ****/
/***                     not interested in (e.g. solvent water) to be reduced   ****/
/***                     to just the solvate for analysis with appropriate      ****/
/***                     alterations to the FIELD file.                         ****/

         if ( natms > num_atoms+1)
           {
             printf("The FIELD file gives a smaller number so will assume the HISTORY file\n");
             printf("has been reduced.\n");
           }
         else
           {
             printf("\nERROR: The FIELD file gives a larger number, so have to exit.\n");
             exit(0);
           }
       }

/**** trajectory key have we co-ordinates, velocities etc ***/
     keytrj    = get_int(&ichar[0], &place, &itsanum,
                         &ndigi, num_of_chars, &sign);

/**** trajectory key have we co-ordinates, velocities etc ***/
     imcon     = get_int(&ichar[0], &place, &itsanum,
                         &ndigi, num_of_chars, &sign);

/**** integration time step  ********************************/
      *p_tstep = get_doub( &ichar[0], num_of_chars, &place, &itsanum );

   }

 if (pbc > 0)
   {
/**** read in cell information ******************************/
      fscanf(file_fp, "%le %le %le", p_latt_vec  , p_latt_vec+1, p_latt_vec+2);
      fscanf(file_fp, "%le %le %le", p_latt_vec+3, p_latt_vec+4, p_latt_vec+5);
      fscanf(file_fp, "%le %le %le", p_latt_vec+6, p_latt_vec+7, p_latt_vec+8);
   }

/*** False read to get to end of line ****/
     num_of_chars= read_line(file_fp, &ichar[0]);
      if (DEBUG)
        {
          printf("READ 2: ");
          for (iloop=0; iloop < num_of_chars; iloop++) putchar(ichar[iloop]); 
          printf("\n");
        }

/**** read in atoms *****************************************/
for (iloop=0; iloop<= num_atoms; iloop++)
  {
     num_of_chars= read_line(file_fp, &ichar[0]); 

     int_to_string(&ichar[0], &cchar[0], num_of_chars, FALSE );
     sscanf(cchar, "%s", p_molecule->label);

/* Strip any underscores off for pot  **/
     for (ic=0; ic <= strlen(cchar); ic++)
                              if (cchar[ic]=='_') cchar[ic]='\0';
     sscanf(cchar, "%s", p_molecule->pot);

/* Look up element symbol */

     if (strncmp(p_molecule->pot, "c",1) == 0 || strncmp(p_molecule->pot, "C",1) == 0 )
       {
         sprintf(p_molecule->elem,"C");
       }
     else if (strncmp(p_molecule->pot, "h",1) == 0 || strncmp(p_molecule->pot, "H",1) == 0)
       {
         sprintf(p_molecule->elem,"H");
       }
     else if (strncmp(p_molecule->pot, "o", 1) == 0 || strncmp(p_molecule->pot, "O", 1) == 0)
       {
         sprintf(p_molecule->elem,"O");
       }
     else if (strncmp(p_molecule->pot, "n", 1) == 0 || strncmp(p_molecule->pot, "N", 1) == 0)
       {
         sprintf(p_molecule->elem,"N");
       }
     else if (strncmp(p_molecule->pot, "az", 2) == 0 || strncmp(p_molecule->pot, "Al", 2) == 0)
       {
         sprintf(p_molecule->elem,"Al");
       }
     else if (strncmp(p_molecule->pot, "si", 2) == 0 || strncmp(p_molecule->pot, "Si", 2) == 0)
       {
         sprintf(p_molecule->elem,"Si");
       }
     else if (strncmp(p_molecule->pot, "s", 1) == 0 || strncmp(p_molecule->pot, "S", 2) == 0)
       {
         sprintf(p_molecule->elem,"S");
       }
     else if (strncmp(p_molecule->pot, "f", 1) == 0 || strncmp(p_molecule->pot, "F", 2) == 0)
       {
         sprintf(p_molecule->elem,"F");
       }
     else 
       {
          printf("Unable to convert label >>%s<< with pot >>%s<< to element in read_hist_frame.\n",
                                     p_molecule->label, p_molecule->pot);
          printf("ERROR: Atomtype not founded\n");
          exit(0);
          
       }

/*** Assign Neutron scattering factors, i.e. call atomic_bscat_list routine ***/
/*** If a number -100.0 is returned print an error and exit                 ***/

p_molecule->bscat = atomic_bscat_list(p_molecule->elem);

if (p_molecule->bscat < -50.0)
  {
     printf("ERROR : No neutron scattering data available for element %s, update data.h information\n", p_molecule->elem);
     exit(0);
  }

/* printf("For element %s assigning bscat %10.6f\n", p_molecule->elem, p_molecule->bscat); */

/*** End of new Neutron scattering factor assignment **************************/

/* falsify group number and name ***/

     sprintf(p_molecule->group,"DLPY");
     sprintf(p_molecule->group_no,"1");

     place= 7;
     dummy= get_int(&ichar[0], &place, &itsanum,
                         &ndigi, num_of_chars, &sign);
  
     if (itsanum)
       {
          p_molecule->mass= get_doub( &ichar[0], num_of_chars, &place, &itsanum ); 
       }
     else
       {
       }

     if (itsanum)
       {
          p_molecule->part_chge= get_doub( &ichar[0], num_of_chars, &place, &itsanum );
       }


/*****************************/
/* Read in atom coordinates  */
/*****************************/

     fscanf(file_fp, "%le %le %le", &(p_molecule->x), &(p_molecule->y), &(p_molecule->z));
     if (levcfg > 0)
       {
          fscanf(file_fp, "%le %le %le", &(p_molecule->vx), &(p_molecule->vy), &(p_molecule->vz));
       }
     if (levcfg > 1)
       {
          fscanf(file_fp, "%le %le %le", &force[0], &force[1], &force[2]);
       }

     p_molecule++;
/*** False read to get to end of line ****/
     num_of_chars= read_line(file_fp, &ichar[0]); 
/*      if (DEBUG)
        {
          printf("READ false: ");
          for (iloop=0; iloop < num_of_chars; iloop++) putchar(ichar[iloop]); 
          printf("\n");
        } */
  }
DEBUG=FALSE;

/***** Need to read another line to get past the \n of the last line *****/
/***** and so prepare for the next frame read ****************************/

/*if (levcfg == 0) num_of_chars= read_line(file_fp, &ichar[0]); */

return 10;
}

