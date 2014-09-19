/*********************************************************/
/***** Read the information about nonbond parameters *****/
/***** from lines begining with @ !                  *****/
/***** Dave Willock  November 1995                   *****/
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"
#include "data.h"

int locate_string_in_string( char *p_key, char *p_char, int num_of_chars );

int string_from_int(int *p_int, char *p_string);

int next_none_space( int *p_ichar, int start, int num_of_chars );

int next_space( int *p_ichar, int start, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

void get_nonbond_params(FILE *fp, int *p_line, char *p_line_string, 
                        int *p_num_of_chars)
{
#include "header.h"
int ipoint, iend, itsanum;
int idummy, ref;
int num_digi, sign;

double version;

/**************************************************************/
/****** Ignore titles for parameters                    *******/
/**************************************************************/

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
printf("Entered with >>%s<< that has %d characters\n",p_line_string, *p_num_of_chars);
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/


   while ( locate_string_in_string(TITLE_LINE,       p_line_string, *p_num_of_chars)
        || locate_string_in_string(ILLUSTRATION_LINE, p_line_string, *p_num_of_chars)
        || locate_string_in_string(JUST_RETURN      , p_line_string, *p_num_of_chars)
        || *p_num_of_chars == 0 || *p_num_of_chars == -1)
     {
        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
printf("Looking at >>%s<<\n",p_line_string);
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
     }

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
printf("Looking at >>%s<<\n",p_line_string);
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/

/**************************************************************/
/****** Stay in the loop till we run out of parameter   *******/
/****** lines                                           *******/
/**************************************************************/

   num_potential_types=-1;
   while(*p_num_of_chars > 0)
     {

        num_potential_types++;
 
        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (potent[num_potential_types].pot, 
                                    p_line_string+ipoint, iend-ipoint);

        ipoint= iend;
        potent[num_potential_types].a
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        potent[num_potential_types].b= 
                   get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
   
         
/***************************************************************/
/****** Read in next parameter line  **************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
return;
}
