#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int find_field(int *p_ichar, int field_desired, int start, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

void read_rdf(FILE *file_fp, int *p_num_molecules, int *p_num_this_mol,
                list_partition *p_demarcation, int *p_num_types )
{

int good_read, num_of_chars, itsanum;
int ichar[LINESIZ], ndigi, sign;

int is_molecules, is_num_this_mol, is_demarcation;
int last_start, place, num, occurances, iloop;
int first_time;

char *p_key;

/*********************************************************************/
/*** Prepare the ground **********************************************/
/*********************************************************************/

last_start= 0;
*p_num_types = -1;
first_time= FALSE;

/*********************************************************************/
/**** Read in the bits we want from the FIELD file *******************/
/*********************************************************************/

good_read= TRUE;

while (good_read)
  {
    num_of_chars= read_line(file_fp, &ichar[0]);

    good_read= num_of_chars != -10;   

    if (good_read)
      {
         p_key= "MOLE";
         is_molecules= locate_string( p_key, &ichar[0], num_of_chars);
  
         p_key= "NUMM";
         is_num_this_mol= locate_string( p_key, &ichar[0], num_of_chars);

         p_key= "ATOM";
         is_demarcation= locate_string( p_key, &ichar[0], num_of_chars);

         place=0;
         if (is_molecules && !first_time )
           {

             first_time = TRUE;
             *p_num_molecules= get_int( &ichar[0], &place, &itsanum, &ndigi, num_of_chars, &sign );
             printf("In read_field found %d molecules\n", *p_num_molecules);

             if (!itsanum)
               {
                 printf("Error: Missing number on MOLECULES line in FIELD file.\n");
               } 
           }
         else if (is_num_this_mol)
           {
             ++*p_num_types;
             occurances = get_int( &ichar[0], &place, &itsanum, &ndigi, num_of_chars, &sign );
             if (!itsanum)
               {
                 printf("Error: Missing number on NUMMOLS line in FIELD file.\n");
                 occurances= 0;
               }
            
             *p_num_this_mol = occurances;
             p_num_this_mol++;

           }
         else if (is_demarcation)
           {

             num = get_int( &ichar[0], &place, &itsanum, &ndigi, num_of_chars, &sign );
             if (!itsanum)
               {
                 printf("Error: Missing number on ATOMS line in FIELD file.\n");
                 num=0;
               }

             for (iloop=0; iloop < occurances; iloop++)
               {
                 p_demarcation->start = last_start;
                 p_demarcation->num   = num;
                 p_demarcation->end   = p_demarcation->start + num - 1;
                 last_start           = p_demarcation->end +1;
                 p_demarcation++;
               }
           }
      }
  }

fclose(file_fp);

return;

}

