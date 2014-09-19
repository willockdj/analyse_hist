#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void do_mopac(FILE *file_fp, atom *p_molecule, int num_atoms, int need_grads )
{

int iloop, dummy;

char call_mopac[FILELEN_MAX];

/*********************************************************************/
/**** Write out MOPAC input file *************************************/
/*********************************************************************/

if (need_grads)
  {
     fprintf(file_fp, "XYZ MMOK PM3 1SCF GRADIENTS MULLIK GEO-OK\ndriven file\ndriven file\n");
  }
else
  {
     fprintf(file_fp, "XYZ MMOK PM3 1SCF MULLIK GEO-OK\ndriven file\ndriven file\n");
  }

for (iloop=0; iloop < num_atoms; iloop++)
  {
    fprintf(file_fp, "%s   %16.13f  0  %16.13f  0  %16.13f  0\n",
                  p_molecule->elem, p_molecule->x, p_molecule->y, p_molecule->z);
    p_molecule++;
  }

fclose(file_fp);
sprintf(call_mopac,"mopac_job\n");
dummy= system(call_mopac);

return;

}

