/***********************************************************************/
/**print_energy.c                                                      */
/*                                                                     */
/* Taken from output_routines in Zebedde, Dave and Dewi April 1997     */
/***********************************************************************/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void print_dashes(int ndashes,FILE *fp);

void print_energy(FILE *fp, energy *p_energy, atom *p_molecule, 
					int num_atoms, int madalung)
{

#include "header.h"
int iatom;
atom *p_atom;

print_dashes(60,fp);
fprintf(fp, "Internal Energies for Template\n");
fprintf(fp,"Stretch energy     : %lf\n",intra_energy.stretch);


fprintf(fp, "Template host interaction energy\n");
fprintf(fp,"Steric Energy Flag : %i\n",p_energy->steric);
fprintf(fp,"Non-bonded Energy  : %lf\n",p_energy->non_bonded);
fprintf(fp,"Coulombic Energy   : %lf (Template->host only= %lf)\n",p_energy->charges, p_energy->guest_host);
fprintf(fp,"Total interaction  : %lf\n", p_energy->non_bonded+p_energy->charges);
fprintf(fp,"\nEnergies from Minimizer (last run)\n");
fprintf(fp,"\t\t\t\t\t\t\tinitial\t\tend\n");
fprintf(fp,"Total Energy      : \t%10.5f\t%10.5f\n",
					p_energy->minimizer_init_total,
					p_energy->minimizer_end_total);
fprintf(fp,"Non-bonded Energy : \t%10.5f\t%10.5f\n",
					p_energy->minimizer_init_nonbond,
					p_energy->minimizer_end_nonbond);
print_dashes(60,fp);

if (madalung)
  {

      fprintf(fp,"Atomic contributions to non-bond potential\n\n");

      fprintf(fp,"Atom No  type  vdw energy\n\n");

      p_atom= p_molecule-1;
      for (iatom=0; iatom <= num_atoms; iatom++)
         {
           p_atom++;
           fprintf(fp,"  %-4d    %-4s    %10.6f\n",
                                 iatom+1, p_atom->label, p_atom->vdw_energy); 
         }

      print_dashes(60,fp);
   }

return;
}
