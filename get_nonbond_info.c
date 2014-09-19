/*********************************************************/
/***** Read the information about nonbond parameters *****/
/***** from lines begining with @ !                  *****/
/***** Dave Willock  November 1995                   *****/
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_values.h"
#include "data.h"

int locate_string_in_string( char *p_key, char *p_char, int num_of_chars );

int next_none_space( int *p_ichar, int start, int num_of_chars );

int next_space( int *p_ichar, int start, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

int string_from_int(int *p_int, char *p_string);

void get_nonbond_info(FILE *fp, 
                      int *p_line, char *p_line_string,
                      int *p_num_of_chars)
{
char grep_string[81],grep_second_string[81];

int ipoint, iend, idummy;

/**************************************************************/
/****** Stay in the loop till we run out of information *******/
/****** lines                                           *******/
/**************************************************************/

   while(locate_string_in_string(INFO_LINE, p_line_string, *p_num_of_chars))
     {
/**************************************************************/
/***** Search the information line for known titles ***********/
/**************************************************************/

       if (locate_string_in_string(POT_COMBINATION, p_line_string, *p_num_of_chars))
         {
           ipoint= next_space(p_line, 0, *p_num_of_chars);
           ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
           iend= next_space(p_line, ipoint, *p_num_of_chars);

           strncpy (pot_info.combination, p_line_string+ipoint, iend-ipoint);
         }


       else if (locate_string_in_string(POT_TYPE, p_line_string, *p_num_of_chars))
         {
           ipoint= next_space(p_line, 0, *p_num_of_chars);
           ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
           iend= next_space(p_line, ipoint, *p_num_of_chars);

           strncpy (pot_info.type, p_line_string+ipoint, iend-ipoint);
         }

/***************************************************************/
/****** Read in next information line **************************/
/***************************************************************/

           *p_num_of_chars=read_line(fp, p_line);
           idummy= string_from_int(p_line, p_line_string);
     }
return;
}
