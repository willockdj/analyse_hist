/*************************************************************/
/* work out the radius of gyration for the molecule **********/
/* assuming centre of mass and total mass already provided ***/
/* Dave Willock, March 2011.                               ***/
/*************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

#define DEBUG FALSE

/* protype list for this routine */

void min_image( double *x, double *y, double *z);

/*---------------------------------------------------------------------------*/

void radius_gyration(atom *p_molecule, int num_atoms, double *p_c_of_m, 
                     double total_mass, double *p_rgyr )
{
  int icomp, iatom;
  atom *p_atom;

  double vec[3], dist2;

  double cubic_coeffs[3];

/****** Zero rgyr value *******/

   *p_rgyr=0.0;

   p_atom=p_molecule;
   for (iatom=0; iatom<=num_atoms; iatom++)
     {
/*** get atom to centre of mass vector ***/
       vec[0] = p_atom->x - *p_c_of_m;
       vec[1] = p_atom->y - *(p_c_of_m+1);
       vec[2] = p_atom->z - *(p_c_of_m+2);

       min_image( &vec[0], &vec[1], &vec[2]);

       dist2= vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];

       *p_rgyr += p_atom->mass * dist2;

       p_atom++;
     }

   *p_rgyr = sqrt(*p_rgyr / total_mass);

   if (DEBUG) printf("Make this molecule's radius of gyration %10.6f Angstroms\n", *p_rgyr);

return;
}
