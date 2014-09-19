/******************************************************************/
/*** make_msd_traj is a routine to ensure that the minimum image **/
/*** trajectory from a HISTORY file is unravelled into a real    **/
/*** trajectory for msd calculation. This is done by working out **/
/*** if the first atom in the list has crossed a periodic        **/
/*** boundary since the last frame, if it has it is moved to the **/
/*** position outside the unit cell it should occupy. All other  **/
/*** atoms in the molecule are then moved to a minimum image     **/
/*** position with respect to the first.                         **/
/*** Note that in this scheme only the first atom requires       **/
/*** image_disp and last_image data .                            **/
/*** Dave Willock, Cardiff August 2005.                          **/
/******************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"
#include "header.h"

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void min_image( double *x, double *y, double *z);

void make_msd_traj(atom *p_molecule, int start_mol, int end_mol,
                   atom *p_last_frame, vector *p_image_disp, int frame_index)
{

int iatom;

atom *p_atom;

double vec[3];
double dx, dy,dz, da, db, dc;

/*******************************************************************/
/** Two cases: first call work out initial displacement of first ***/
/**            atom.                                             ***/
/**            second call update first atom position,           ***/
/*******************************************************************/
if (frame_index == 0)
 {
    printf("Processing first frame, start_mol=%d end_mol = %d\n", start_mol, end_mol);

/* Set start point as zero image displacement  */
                                           
    p_image_disp->x = 0.0;
    p_image_disp->y = 0.0;
    p_image_disp->z = 0.0;
 }
else 
 {

/***********************************************************/
/* Work out displacement of first atom from last time step */
/* Effective current position of atom is current co-ords   */
/* plus image_disp vector                                  */
/***********************************************************/

    dx = p_last_frame->x - p_molecule->x - p_image_disp->x;
    dy = p_last_frame->y - p_molecule->y - p_image_disp->y;
    dz = p_last_frame->z - p_molecule->z - p_image_disp->z;

    cart_to_fract( dx, dy, dz, &da, &db, &dc, &recip_latt_vec[0] );

/* Check for lattice vector size jumps that indicate change of min-image */

    if( da > 0.5)
      {
         p_image_disp->x += latt_vec[0];
         p_image_disp->y += latt_vec[1];
         p_image_disp->z += latt_vec[2];
      }
    else if( da <-0.5)
      {
         p_image_disp->x -= latt_vec[0];
         p_image_disp->y -= latt_vec[1];
         p_image_disp->z -= latt_vec[2];
      }

    if( db > 0.5)
      {
        p_image_disp->x += latt_vec[3];
        p_image_disp->y += latt_vec[4];
        p_image_disp->z += latt_vec[5];
      }
    else if( db < -0.5)
      {
        p_image_disp->x -= latt_vec[3];
        p_image_disp->y -= latt_vec[4];
        p_image_disp->z -= latt_vec[5];
      }

    if( dc > 0.5)
      {
        p_image_disp->x += latt_vec[6];
        p_image_disp->y += latt_vec[7];
        p_image_disp->z += latt_vec[8];
      }
    else if( dc < -0.5)
      {
        p_image_disp->x -= latt_vec[6];
        p_image_disp->y -= latt_vec[7];
        p_image_disp->z -= latt_vec[8];
      }
    }

    p_molecule->x += p_image_disp->x;
    p_molecule->y += p_image_disp->y;
    p_molecule->z += p_image_disp->z;

    *p_last_frame = *p_molecule;

    p_atom = p_molecule;
    for (iatom=start_mol+1; iatom <= end_mol; iatom++)
      {
/* Get vector from current atom to first in list */

         p_atom++;
         vec[0]= p_atom->x - p_molecule->x;
         vec[1]= p_atom->y - p_molecule->y;
         vec[2]= p_atom->z - p_molecule->z;

        /* Find minimum image of this atom from first in list */

         min_image(&vec[0], &vec[1], &vec[2]); 

        /* Move the atom to its minimum image co-ordinates */

         p_atom->x = p_molecule->x+vec[0];
         p_atom->y = p_molecule->y+vec[1];
         p_atom->z = p_molecule->z+vec[2];
      }
return;
}
