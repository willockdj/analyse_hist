/***************************************************************/
/* gather_molecule.c : move atoms within a molecule so that it */
/*                     is correctly assembled with no pbc      */
/*                     breaks.                                 */
/*                     routine added to avoid problems of      */
/*                     molecules that are large compared to    */
/*                     half the periodic cell dimension.       */
/*                     Dave Willock July 09.                   */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B);

void min_image( double *x, double *y, double *z);

void gather_molecule(atom *p_molecule, int num_atoms, int which_mol)
{
#include "header.h"
int ineigh, iflag, num_ends, atoms_at_ends[MAX_ENDS];
int flag[MAX_ATOMS];
int this_neigh, iatom, iends;
double vec[3];
atom *p_current_atom, *p_atom;
atom *p_atom1, *p_atom2, *p_neigh;

/********* zero the flags ************/

for (iflag=0; iflag<= num_atoms; iflag++) flag[iflag]=FALSE;

/*** Start search with atom zero for this molecule ***/

num_ends= 0;

if (which_mol < 0)
  {
    atoms_at_ends[0]= 0;
    flag[0]=TRUE;
  }
else
  {
    p_atom = p_molecule;
    for (iatom=0; iatom<= num_atoms; iatom++) 
       {
         if ( p_atom->mol == which_mol)
           {
              flag[iatom]=TRUE;
              break;
           }
         p_atom++;
       }

    if (iatom == num_atoms+1)
      {
         printf("ERROR: cannot find molecule %d in gather_molecule\n", which_mol);
         exit(0);
      }

    atoms_at_ends[0]= iatom;
    flag[iatom]=TRUE;
  }

/******* set the flags for all end atoms and generate new ends  ******/
/******** till we have been everywhere in the molecule          ******/

while (num_ends >= 0)
  {
    p_current_atom= p_molecule + atoms_at_ends[0];

    for (ineigh = 0; ineigh <= (p_current_atom->num_neigh); ineigh++)
      {
        this_neigh=  p_current_atom->neighb[ineigh];

/***************************************************************************/
/*** Shift the neighbours to min_image with respect to this atom ***********/
/***************************************************************************/

        p_neigh= p_molecule+this_neigh;

        join_atoms( p_current_atom, p_neigh, &vec[0]);

/*** take min_image and shift atom2 to it ****************************/

        min_image(&vec[0], &vec[1], &vec[2]);

        p_neigh->x = p_current_atom->x + vec[0];
        p_neigh->y = p_current_atom->y + vec[1];
        p_neigh->z = p_current_atom->z + vec[2];

/***************************************************************************/
/***** add this neighbour to the ends list if it is not already flagged ****/
/***** dissallow atom1 joining the list to stop the other chunk being   ****/
/***** crossed into                                                     ****/
/***************************************************************************/

        if ( !flag[this_neigh] )
          { 
            flag[this_neigh]= TRUE;
            if (num_ends < MAX_ENDS && (p_molecule+this_neigh)->num_neigh != 0) 
              {
                num_ends++;
                atoms_at_ends[num_ends]= this_neigh;       
              }
            else if (num_ends > MAX_ENDS)
              {
                printf("Run out of ends in find_chunk.c increase MAX_ENDS");
                exit(1);
              }
          }
      } 

/******* shuffle atoms_at_ends list to cover the one we have dealt with ****/

    for (iends= 0; iends < num_ends; iends++)
      {
          atoms_at_ends[iends] =  atoms_at_ends[iends+1];   
      }
    num_ends--;
  }

return;
}

		
