/***************************************************************************/
/*** Apply a true square matrix operator to a set of atoms *****************/
/*** started Dave Willock May 10th 1997                    *****************/
/***************************************************************************/
#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void apply_matrix(atom *p_molecule, double *p_matrix, int num_atoms)
{

  int iloop,icoloum;

  double hold_coord[3], *p_elem;

  atom *p_atom;

printf("DEBUG>> Have %d atoms in apply_matrix\n",num_atoms);

  for (iloop= 0; iloop <= num_atoms; iloop++)
  {
     p_elem= p_matrix;
     p_atom= p_molecule+iloop; 
     for (icoloum=0; icoloum< 3; icoloum++)
     {
        hold_coord[icoloum]= 0.0;
        hold_coord[icoloum]+= *p_elem*(p_atom->x);
        p_elem++;
        hold_coord[icoloum]+= *p_elem*(p_atom->y);
        p_elem++;
        hold_coord[icoloum]+= *p_elem*(p_atom->z);
        p_elem++;
     } 
     p_atom->x = hold_coord[0];
     p_atom->y = hold_coord[1];
     p_atom->z = hold_coord[2];
  }

 return;
}

