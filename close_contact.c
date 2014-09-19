#include <stdio.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"
#include "global_values.h"

double atom_separation_squared(atom *p_A, atom *p_B, int pbc);

void close_contact(atom *p_molecule1, int num_atoms1,
                   atom *p_molecule2, int num_atoms2,
                   double cutoff2, int pbc, int *p_num_found,
                   pair_list *p_clash_list )
{
int iloop, jloop;

atom *p_atom1, *p_atom2;

double dist;

*p_num_found= -1;
p_atom1= p_molecule1;
for (iloop=0; iloop <= num_atoms1; iloop++)
  {
    p_atom2= p_molecule2;
    for (jloop=0; jloop <= num_atoms2; jloop++)
      {
        if (iloop != jloop)
          {
            dist = atom_separation_squared(p_atom1, p_atom2, pbc); 

            if (dist <= cutoff2)
              {
                ++*p_num_found;

                strcpy(&(p_clash_list->label1[0]), &(p_atom1->label[0])); 
                strcpy(&(p_clash_list->label2[0]), &(p_atom2->label[0])); 

                p_clash_list->index1= iloop;
                p_clash_list->index2= jloop;

                p_clash_list->sep= sqrt(dist);

                p_clash_list++;
              } 
            p_atom2++;
          }
       }
    p_atom1++;
  }

return;
}
