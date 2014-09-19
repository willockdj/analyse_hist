/****************************************************************************************/
/***** Simple summation electrostatics for clusters and removing self interaction *******/
/***** from interaction energy calcs ****************************************************/
/***** Since we only care about the basic molecule only molecule A is assigned an *******/
/***** energy                                                                     *******/
/***** Started 25 July 1996 Dave Willock ************************************************/
/****************************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "structures.h"
#include "ewald.h"
#include "maxima.h"
#include "own_maths.h"

void isol_coul(atom *p_mol_A, int num_mol_A_atoms, 
               atom *p_mol_B, int num_mol_B_atoms, int zero_first, int self)
{
#include "header.h"
int mol_index_A, mol_index_B;
atom *p_atom_A, *p_atom_B;
double dx, dy, dz, separation;

if (zero_first) 
  {
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
         p_atom_A= p_mol_A+mol_index_A;
         p_atom_A->electrostatic_pot = 0;
      }
  }


/**************************************************************************************/
/****** work out space sum limits needs implementing for cluster work!!! **************/
/**************************************************************************************/

/* num_needed[0] = coul_sum_max/real_latt_sizes[0]; */
/* num_needed[1] = coul_sum_max/real_latt_sizes[1]; */
/* num_needed[2] = coul_sum_max/real_latt_sizes[2]; */

/***************************************************************************************/
/**** Case when mol_A not mol_B ********************************************************/
/***************************************************************************************/

if (!self)
  {
    printf ("Cluster electrostatic sum not yet available!!!\n");
    exit(0);
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
        p_atom_A= p_mol_A+mol_index_A;
        for (mol_index_B=0; mol_index_B <= num_mol_B_atoms; mol_index_B++)
          {
            p_atom_B= p_mol_B+mol_index_B;
          }
      }
  }

/***************************************************************************************/
/**** Case when we wish to remove mol_A/mol_A interactions from the sum ****************/
/***************************************************************************************/

else
  {
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
        p_atom_A= p_mol_A+mol_index_A;
        for (mol_index_B=mol_index_A+1; mol_index_B <= num_mol_B_atoms; mol_index_B++)
          {
            p_atom_B= p_mol_B+mol_index_B;

            dx= (p_atom_A->x)-(p_atom_B->x);
            dy= (p_atom_A->y)-(p_atom_B->y);
            dz= (p_atom_A->z)-(p_atom_B->z);
  
            separation = sqrt (dx*dx + dy*dy + dz*dz);

            p_atom_A->electrostatic_pot -= (p_atom_B->part_chge)/separation;
          }
      }
  }

return;
}
