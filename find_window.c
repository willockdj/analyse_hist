/***************************************************************/
/* find_window.c: look for the atoms defining a window in a   **/
/******           molecule.                                   **/
/******     started April 2011 Dave Willock                   **/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

double size_vector(double *p_vector);

void min_image( double *x, double *y, double *z);

int neighbour_order(atom *p_molecule, int num_atoms, int atom1, int atom2);

int find_window(atom *p_molecule, int num_atoms, windef *p_winref, winlist *p_winsets )
{
#include "header.h"
int iii, jjj, ineigh, iflag, iwin;
int this_neigh, iatom, iends, order;
int num_matched, have_matched, found;
int matched_indices[MAX_WINATOMS_PER_MOL];
int matched_neighs[MAX_WINATOMS_PER_MOL];
int got[MAX_WINATOMS_PER_MOL];
atom *p_neigh, *p_atom;

double vec1[3], vec2[3], dot;

found=FALSE;
printf("Looking for a window or two.....\n");

/**** Loop over molecule looking for atoms that match definition of window atoms ***/

num_matched=-1;
p_atom=p_molecule;
for (iii=0; iii<=num_atoms; iii++)
   {
     have_matched=FALSE;

/*** Look over atoms in reference definition ****/

     for (iwin= 0; iwin < p_winref->num_atoms; iwin++) 
       {
         if (!have_matched && strcmp(p_atom->label, p_winref->atoms[iwin]) == 0)
           {
             if (DEBUG) printf("Atom %d has the right label %s (%s) it has %d neighbours to check\n", iwin,
                                           p_atom->label, p_winref->atoms[iwin], p_atom->num_neigh); 

/*** Look over the neighbours of the label matched atom ****/

             for (ineigh=0; ineigh<= p_atom->num_neigh; ineigh++)
                {
                   p_neigh=p_molecule+p_atom->neighb[ineigh];

                    if (strcmp(p_neigh->label, p_winref->neigh[iwin]) == 0)
                      {
/**** DEBUG CHECK START **/
                         vec2[0] = p_atom->x - p_neigh->x;
                         vec2[1] = p_atom->y - p_neigh->y;
                         vec2[2] = p_atom->z - p_neigh->z;

                         min_image( &vec2[0], &vec2[1], &vec2[2]);
/**** DEBUG CHECK END ****/
                         if (DEBUG) printf(".....it also has the right neighbour %s (%s) they are %10.6f A apart\n",
                                           p_neigh->label, p_winref->neigh[iwin], size_vector(&vec2[0])); 
/******************************************************************/
/*** Found candidate add its index to the matched_indicies list ***/
/*** and set the corresponding got[] flag FALSE.                ***/
/******************************************************************/
                         have_matched=TRUE;
                         num_matched++;
   
                         if (num_matched < MAX_WINATOMS_PER_MOL)
                           {
                             matched_indices[num_matched]=iii; 
                             matched_neighs[num_matched]=p_atom->neighb[ineigh]; 
                             got[num_matched]=FALSE;
                             found=TRUE;
                           }
                         else
                           {
                             printf("ERROR: Found too many window defining atoms, maximum allowed : %d\n",
                                                                MAX_WINATOMS_PER_MOL);
                             exit(0);
                           }
                      }
                }
           }
       }

     p_atom++;
   }

printf("Matched %d atoms as possible window atoms\n", num_matched+1);

/*** Find sets of window atoms for each window in the molecule                      ***/
/*** This depends on the neighbour order and the dot products for neighbour -> atom  **/
/*** bond. Currently this is hard coded for Cage 3 so that neighbour order has to be **/
/*** 11 and the bond vectors have to be aligned with a positive dot product.         **/

/* p_winsets->num_windows=-1; */

/*** loop over matched atoms ****/
for (iii=0; iii<=num_matched; iii++)
   {
      if (!got[iii])
        {
              printf("DEBUG >> not got...\n");
          (p_winsets->num_windows)++;
          iwin= p_winsets->num_windows;

          p_winsets->num_atoms[iwin]=0;

          p_winsets->iatom[iwin][0]=matched_indices[iii];

          got[iii]=TRUE;
              printf("DEBUG >> now got...\n");
          p_atom= p_molecule+matched_indices[iii];
          p_neigh= p_molecule+matched_neighs[iii];
          if (DEBUG) printf("Atom %d matches, label %s as it has neigh %s\n", matched_indices[iii], 
                                                                                  p_atom->label,
                                                                                  p_neigh->label);
          vec1[0] = p_atom->x - p_neigh->x;
          vec1[1] = p_atom->y - p_neigh->y;
          vec1[2] = p_atom->z - p_neigh->z;
          min_image( &vec1[0], &vec1[1], &vec1[2]);

/*** loop over rest of atoms that have been matched looking for ones in the same window ***/
              printf("DEBUG >> looking for partners...\n");
          for (jjj=iii+1; jjj<=num_matched; jjj++)
             {
              if (!got[jjj]) 
               {   
                 printf("DEBUG Not got this potential partner into a window yet....\n");     
                p_atom= p_molecule+matched_indices[jjj];
                p_neigh= p_molecule+matched_neighs[jjj];

                vec2[0] = p_atom->x - p_neigh->x;
                vec2[1] = p_atom->y - p_neigh->y;
                vec2[2] = p_atom->z - p_neigh->z;
                min_image( &vec2[0], &vec2[1], &vec2[2]);
         
                dot = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
                dot = dot / size_vector(&vec1[0]);
                dot = dot / size_vector(&vec2[0]);

                order=neighbour_order(p_molecule, num_atoms, 
                                      matched_indices[iii], matched_indices[jjj]);

                DEBUG=TRUE;
                if (DEBUG)  
                   printf("iii: %d (atom %d)  jjj: %d (atom %d) have dot %10.6f they are %d order neighbours.",
                                                                         iii, matched_indices[iii],
                                                                         jjj, matched_indices[jjj],
                                                                         dot, order);

/****	Altered for cp definition of windows, Mar 12 DJW ***/
/*                if (dot > 0.1 && order == 11)   */
                if (dot > 0.1 && order == 9)  
                  { 
                    got[jjj]=TRUE;
                    if (DEBUG) printf("...... so have a match would use in window %d\n", iwin);

/**** Increment the number of atoms defining this window ****/
/**** Pick up the atom index that is needed              ****/

                   (p_winsets->num_atoms[iwin])++;
                   iatom= p_winsets->num_atoms[iwin];

                   p_winsets->iatom[iwin][iatom]=matched_indices[jjj];
                 }
               else
                 {
                   if (DEBUG) printf("match failed\n");
                 }
                DEBUG=FALSE;
              }
            }
          if (DEBUG) printf("\n");
       }
   }

return found;
}

		
