#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"

/* protype list for this routine */
void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

void gather_molecule(atom *p_molecule, int num_atoms, int which_mol);

double size_vector(double *p_vector);

/*---------------------------------------------------------------------------*/

/* work out furthest atom in the direction of a given vector. */
/* starting at a centre provided.                             */

int furthest_atom(atom *p_molecule, int num_atoms, double *p_centre, double *p_vec, double *p_dist)

{
  int imax, this_atom;
  atom *p_atom;

  double dot, size_v;
  double vec[3];

  double *p_vec1, *p_vec2;

/*----------------------------------------------------------------------------*/

/****** move molecule to centre of mass *****/
  printf("Moving molecule to centre at %10.6f %10.6f %10.6f\n", 
                                                *p_centre, *(p_centre+1), *(p_centre+2)); 
  vec[0]= - *p_centre;
  vec[1]= - *(p_centre+1);
  vec[2]= - *(p_centre+2);
  move_molecule(p_molecule, num_atoms, &vec[0]);

  gather_molecule(p_molecule, num_atoms, -1);

/****** find furthest atom along vector direction ****/

   p_atom= p_molecule;
   p_vec1 = p_vec+1;
   p_vec2 = p_vec+2;
   *p_dist = 0.0;
   imax = -1;
   size_v= size_vector(p_vec);
   
   for (this_atom=0; this_atom <= num_atoms; this_atom++)
     {

        dot =   p_atom->x * *p_vec 
              + p_atom->y * *p_vec1
              + p_atom->z * *p_vec2;

        if (fabs(dot) > *p_dist) 
          {
            *p_dist= fabs(dot);
            imax = this_atom;
          }
       p_atom++;
     }

   if (imax < 1)
     {
       printf("ERROR: Cannot find any suitable atoms in furthest atom rountine....\n");
       exit(0);
     }
/****** move molecule back from centre of mass *****/

   move_molecule(p_molecule, num_atoms, p_centre); 

   *p_dist = *p_dist/size_v;

return imax;
}
