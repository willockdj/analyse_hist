/*********************************************************/
/***** Calculate the molecules dipole moment         *****/
/***** Dave Willock  April 1999                      *****/
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_values.h"
#include "maxima.h"
#include "structures.h"
#include "reader.h"
#include "header.h"

int mol_dipole( atom *p_molecule, int num_atoms, vector *p_mol_dipole )
{
atom *p_atom;

double tot_chge;

int error=FALSE;
int iatom;

/**** Calculate charge contributions to dipole ****/

tot_chge=0.0;
p_atom=p_molecule;
p_mol_dipole->x = 0.0;
p_mol_dipole->y = 0.0;
p_mol_dipole->z = 0.0;

for (iatom=0; iatom <= num_atoms; iatom++)
   {
      tot_chge+= p_atom->part_chge;

      p_mol_dipole->x += p_atom->x * p_atom->part_chge;
      p_mol_dipole->y += p_atom->y * p_atom->part_chge;
      p_mol_dipole->z += p_atom->z * p_atom->part_chge;
      p_atom++;
   }

if ( tot_chge > 1E-4 || tot_chge < -1E-4 ) 
  {
    printf("ERROR: Attempt to get dipole moment for charged molecule\n");
    printf("Tot charge: %10.6f\n\n", tot_chge);
    exit(0);
  }

return error;
}
