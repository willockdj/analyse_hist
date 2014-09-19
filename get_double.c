#include <stdio.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

double get_double(int *p_error, int skip_lines)
{
#include "header.h"
double answer;
char *tok;
 
  tok= tok_get(input_fp, skip_lines, FALSE);

/*******************************************************************/
/**** Check all is well, tok will give FALSE if it is NULL *********/
/**** If there is a string there check the first character *********/
/**** to see if it is a valid decimal number               *********/
/**** Dave Willock March 1997 **************************************/
/*******************************************************************/

  if (tok)
    {
       if ((*tok >= '0' && *tok <= '9') || *tok == '-' || *tok == '.')
          {
             answer = atof(tok);
             *p_error = FALSE;
          }
       else
          {
             answer= 0;
             *p_error = TRUE;
          }
    }
  else
    {
       answer= 0;
       *p_error = TRUE;
    }
  return (answer);
}

