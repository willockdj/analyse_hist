#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

void read_statis(FILE *file_fp, double *p_enthalpy, 
                 double *p_time, int *p_num )
{
int num_data;
int iloop;
int done;

char *p_key;
char dummy[100];

/*********************************************************************/
/*** Prepare the ground **********************************************/
/*********************************************************************/

*p_num=-1;

/*********************************************************************/
/**** Read in the bits we want from the STATIS file title ************/
/*********************************************************************/

done = FALSE;
while ( strncmp(dummy,"ENERGY",5) != 0 && !done )
                             done = fscanf( file_fp, "%s",dummy) == EOF;

printf("After energy dummy in STATIS is %s\n",dummy);

if ( done )
  {
     printf("ERROR: No ENERGY in STATIS file\n");
     exit(0);
  }

/**** Skip Energy unit lines ****/
for ( iloop = 0; iloop < 4; iloop++ )
                             done = fscanf( file_fp, "%s",dummy) == EOF;

printf("At start of while dummy is %s\n",dummy);
while ( !done )
  {

/*** skip nstep ***/
    done = fscanf( file_fp, "%s", dummy) == EOF;

    if (!done)
     {

       (*p_num)++;

/*** Read in the time ***************/
       fscanf( file_fp, "%le", p_time);

/*       if (*p_num == 0)   */
/*         {   */
/*            *p_time = atof(dummy); */
/*         } */
/*       else */
/*         { */
/*         } */

       printf("dummy=%s, Statis time %10.6f\n", dummy, *p_time);
       p_time++;

/*** Read in number of following fields *****/

       fscanf( file_fp, "%d", &num_data);

    for ( iloop = 0; iloop < num_data; iloop++ )
     {
       if (iloop==9)
         {
           fscanf( file_fp, "%le", p_enthalpy );
           printf("Statis enthalpy %10.6f\n", *p_enthalpy);
           p_enthalpy++;
         }
       else
         {
           done = fscanf( file_fp, "%s",dummy) == EOF;
         }
      }
    }
  }

return;

}

