#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"
#include "math.h"

void write_distrib_csv(FILE *fp, char *p_title_x, char *p_title_y, 
                       char *p_title_z,
                       double *p_x, int *p_y, int *p_z,
                       int have_z, int num)
  {
int i;

fprintf(fp,"%s,%s",p_title_x,p_title_y); 

if ( have_z )
        fprintf(fp,",%s\n",p_title_z); 

fprintf(fp,"\n"); 


for (i=0; i<=num; i++)
  {
    if (fabs(*p_x) < 1.0e-4 || fabs(*p_x) >1.0e4 )
      {
        fprintf(fp,"%e",*p_x); 
      }
    else
      {
        fprintf(fp,"%f",*p_x); 
      }
    fprintf(fp,","); 

    fprintf(fp,"%d", *p_y); 

    if ( have_z )
      {
        fprintf(fp,","); 

        fprintf(fp,"%d", *p_z); 
        p_z++;
      }
    fprintf(fp,"\n"); 
    p_x++;
    p_y++;
  }

return;
  }

