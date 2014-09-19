/**********************************************************************/
/* generate_neighbours.c : search molecule for neighbour list         */
/* NOTE: neighb is referenced from ZERO                               */
/* started Dave and Dewi 23/3/95                                      */
/**********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"
#include "maxima.h"

#define BOND_TOL 0.1

double standard_bond( char *atom1, char *atom2 );

double atom_separation_squared(atom *p_A, atom *p_B, int pbc);

void generate_neighbours( atom *p_molecule, int num_atoms, 
                          atom_number *p_types, int *p_num_types,
                          int use_pbc, int first_call)
{
#include "header.h"
int iatom1,iatom2;
int idummy2, idummy1;
int itype;

double standard, actual_dist;

atom *p_atom1, *p_atom2;

atom_number *p_this_type;

/********* Work out how many of each atom are present in this molecule *****/

if (first_call)
  {
    *p_num_types= 0;
    strcpy( &(p_types->atom_type[0]), &(p_molecule->elem[0]));
  }

iatom1=0;
p_atom1= p_molecule;

for (iatom1 = 1; iatom1 <= num_atoms; iatom1++)
   {

     p_atom1++;
     p_this_type= p_types-1;

     for (itype=0; itype<= *p_num_types; itype++)
       {
         p_this_type++;
         if (!strcmp(&(p_this_type->atom_type[0]), &(p_atom1->elem[0]))) break;
       }

/**** test to see if that was a known element*****/

     if (itype <= *p_num_types)
       {
         (p_this_type->num)++;
       }
     else
       {
         (*p_num_types)++;
         p_this_type++;
         strcpy( &(p_this_type->atom_type[0]), &(p_atom1->elem[0]));
         p_this_type->num= 0;
       }

   }

/**********************************************************************/
/******************* loop over pairs of atoms in molecule *************/
/*******************         looking for bonds            *************/
/******************* Altered Aug 98 to zero num_neigh for *************/
/******************* all atoms first to allow re-calc of  *************/
/******************* neighbours for an existing molecule  *************/
/******************* Also made num_neigh a bonafide c index ***********/
/******************* i.e. element 0 of neigh list contains  ***********/
/******************* the first neigh index and loops over   ***********/
/******************* neighbours should run <=               ***********/
/******************* Dave Willock                         *************/
/**********************************************************************/

p_atom1= p_molecule;
for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
    p_atom1->num_neigh=-1;
    p_atom1++;
  }

p_atom1= p_molecule-1;
for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
  p_atom1++;
  p_atom2= p_molecule+iatom1;

  for (iatom2=iatom1+1; iatom2 <= num_atoms; iatom2++)
     {
      p_atom2++;

      standard= standard_bond( &(p_atom1->elem[0]), &(p_atom2->elem[0]));

      if (standard != -1) 

        {
           actual_dist = atom_separation_squared(p_atom1, p_atom2, use_pbc);
           actual_dist = sqrt(actual_dist);

/* test bond length against standard */

	   if (actual_dist <= standard+BOND_TOL)
             {
                 (p_atom1->num_neigh)++;
		 idummy1= p_atom1->num_neigh;
	         p_atom1->neighb[idummy1]= iatom2; 

                 (p_atom2->num_neigh)++;
		 idummy2= p_atom2->num_neigh;
	         p_atom2->neighb[idummy2]= iatom1; 
             }

/* end of if (Standard) */
	 }
     }
  } 
  return;
}
