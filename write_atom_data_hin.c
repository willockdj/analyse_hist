
/* routine to write out an atom data line in a hyperchem .hin format */
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"


void write_atom_data_hin(FILE *p_file, atom *p_atom, int index, int start)
{
int iloop, num_neighbours, neigh_index;
int allow[10];

double dist, dx, dy, dz;

/***** Make sure that the atoms do not require pbc to be bonded ***/

  num_neighbours=0; 
  for ( iloop=0; iloop <= p_atom->num_neigh; iloop++)
    {
      neigh_index = p_atom->neighb[iloop] - index+start; 

      dx = p_atom->x - (p_atom+neigh_index)->x;
      dy = p_atom->y - (p_atom+neigh_index)->y;
      dz = p_atom->z - (p_atom+neigh_index)->z;

      dist = sqrt( dx*dx + dy*dy + dz * dz);

      if (dist < 3.0)
        {
          allow[iloop]=TRUE;
          num_neighbours++;
        }
      else
        {
          allow[iloop]=FALSE;
        }

    }

  fprintf( p_file, "atom %d - %2s  ** - 0 %10.7f %10.7f %10.7f %d " ,
           index+1-start,
           p_atom->elem, 
           p_atom->x, p_atom->y, p_atom->z,
           num_neighbours);

  for ( iloop=0; iloop <= p_atom->num_neigh; iloop++)
    {
       if (allow[iloop])
                   fprintf(p_file,"%d s ", p_atom->neighb[iloop]+1);
    }

fprintf(p_file,"\n");
    
return;
}

