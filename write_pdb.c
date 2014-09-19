#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "constants.h"
#include "structures.h"

void write_pdb(FILE *pdb_fp, atom *p_molecule, double *p_abc, int num_atoms)
{
int iloop, iatom, ineigh, neigh_index;
char elem[3];
atom *p_atom;
double dx, dy, dz, dist;
               
/*** VMD requires the correct format of pdb CRYST line as defined at                ****/
/*** http://www.wwpdb.org/documentation/format32/sect8.html                         ****/
/*** Dave Willock July 09                                                           ****/

 fprintf(pdb_fp,"%-8s %6.3f   %6.3f   %6.3f %6.2f %6.2f %6.2f P 1           1\n",
                                            "CRYST1", 
                                            *p_abc, 
                                            *(p_abc+1),
                                            *(p_abc+2),
                                            *(p_abc+3)*RAD_TO_DEG,               
                                            *(p_abc+4)*RAD_TO_DEG,               
                                            *(p_abc+5)*RAD_TO_DEG);              


p_atom = p_molecule;
for ( iloop=0; iloop < num_atoms; iloop++)
  {
    fprintf(pdb_fp,"HETATM %4d  %2s        %4d   %6.3f  %6.3f  %6.3f\n",
                                        iloop+1, p_atom->elem, iloop+1,
                                        p_atom->x, p_atom->y, p_atom->z);
    p_atom++;
  }

p_atom = p_molecule;
for ( iloop=0; iloop < num_atoms; iloop++)
  {
    fprintf(pdb_fp,"CONECT %4d", iloop+1);
     
    for ( ineigh=0; ineigh <= p_atom->num_neigh; ineigh++)
      {
        neigh_index = p_atom->neighb[ineigh]; 

        dx = p_atom->x - (p_molecule+neigh_index)->x;
        dy = p_atom->y - (p_molecule+neigh_index)->y;
        dz = p_atom->z - (p_molecule+neigh_index)->z;

        dist = sqrt( dx*dx + dy*dy + dz * dz);

        if (dist < 3.0) fprintf(pdb_fp," %4d", neigh_index+1);
      }

    fprintf(pdb_fp,"\n");
    p_atom++;
  }

fprintf(pdb_fp,"END\n");

return;
}

