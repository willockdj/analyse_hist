#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

int read_line(FILE *fp, int *p_ichar);

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position, int to_space );

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

int  read_hist_info(FILE *file_fp, char *p_title, int *p_pbc, int *p_levcfg,
                    int have_config, int *p_traj_num ) 
{

int good_read, num_of_chars, itsanum;
int ichar[LINESIZ], ndigi, sign, to_space;

int is_molecules, is_num_this_mol, is_demarcation;
int last_start, place, num, occurances, iloop;
int dummy;

char *p_key;

/*********************************************************************/
/**** Read in the bits we want from the HISTORY file *****************/
/*********************************************************************/

 num_of_chars= read_line(file_fp, &ichar[0]);

 good_read= num_of_chars != -10;   

 if (good_read)
   {
     to_space= FALSE;
     int_to_string(&ichar[0], p_title, num_of_chars, to_space);
/************************************/
/**** Catch case of missing title ***/
/************************************/

     if (strncmp(p_title,"timestep",8) == 0)
       {
         place=0;
         dummy= get_int(&ichar[0], &place, &itsanum,
                             &ndigi, num_of_chars, &sign);

         *p_traj_num= get_int(&ichar[0], &place, &itsanum,
                             &ndigi, num_of_chars, &sign);

         *p_levcfg= get_int(&ichar[0], &place, &itsanum,
                             &ndigi, num_of_chars, &sign);

          *p_pbc   = get_int(&ichar[0], &place, &itsanum,
                             &ndigi, num_of_chars, &sign);

          return 1;
       } 
   }
 else
   {
     printf("Error: No title line included in the HISTORY file\n");
     exit(0);
   }

 num_of_chars= read_line(file_fp, &ichar[0]);

 good_read= num_of_chars != -10;     

 if (good_read)
   {
     place=0;
     *p_levcfg= get_int(&ichar[0], &place, &itsanum,
                             &ndigi, num_of_chars, &sign);

     *p_pbc   = get_int(&ichar[0], &place, &itsanum,
                             &ndigi, num_of_chars, &sign);

/*** Get traj number from timestep lines instead! *****/
/*     if (!have_config) */
/*                 *p_traj_num = get_int(&ichar[0], &place, &itsanum, */
/*                                         &ndigi, num_of_chars, &sign); */
   }
 else
   {
     printf("Error: No periodic boundary information in HISTORY file\n");
     exit(0);
   }


return 0;
}

