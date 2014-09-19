#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"
#include "math.h"

/**************************************************************************/
/** Routine to write out windows information as csv file.               ***/
/** Expects to be called for each type of molecule at each frame        ***/
/**************************************************************************/

void write_window_csv(FILE *fp, double *p_time_now, win_geometry *p_win_rads,
                      int start_type, int num_mols, int num_wins)
  {
int iwin, imol, iwin_tot;
double data;

     iwin_tot=-1;

     printf("In write_window_csv have %d molecules %d windows\n", num_mols, num_wins);
     if (start_type) fprintf(fp,"%f,",*p_time_now); 

     for (imol=0; imol< num_mols; imol++)
        {
           for (iwin=0; iwin< num_wins; iwin++)
              {
                 iwin_tot++;
                 data= p_win_rads->radius[iwin_tot];

                 if (fabs(data) < 1.0e-4 || fabs(data) >1.0e4 )
                   {
                     fprintf(fp,"%e",data); 
                   }
                 else
                   {
                     fprintf(fp,"%f",data); 
                   }

                 if (imol < num_mols-1 || iwin < num_wins-1) fprintf(fp,","); 
              }
        }

return;
  }

