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
/** Routine to write out windows title information to start of  csv file.**/
/** Dave Willock July 2011.                                              **/
/** Debugging June 2013: num_mols here is the number of types of molecules*/
/** in the dlpoly sense and so also have included the use_for_wins array. */
/** Also the num_wins is the total number of windows found for this type  */
/** of molecules so need num_win_per_mol from the winref array in main.   */
/**************************************************************************/

void write_window_csv_titles(FILE *fp,
                             int num_mols, int num_wins, int num_win_per_mol, int *p_use)
  {
int iwin, imol;
int this_mol, this_win;
double data;

fprintf(fp,"time"); 

printf("Arrived in write_window_csv_titles with %d molecules and %d windows, %d per molecule\n",
                    num_mols, num_wins, num_win_per_mol);

     for (imol=0; imol<num_mols; imol++)
        {
           if (*p_use)
              {
                this_mol=1; this_win=1;
                for (iwin=0; iwin<=num_wins; iwin++)
                 {
                   fprintf(fp,","); 

                   fprintf(fp,"m%d_w%d",this_mol, this_win); 

                   this_win++;
                   if (this_win > num_win_per_mol)
                     {
                       this_mol++; this_win=1;
                     }
                 }
             }
          p_use++;
        }
     fprintf(fp,"\n");

return;
  }

