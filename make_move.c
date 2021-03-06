/********************************************************************************/
/**** Routine to make the move suggested            *****************************/
/**** Dave Willock April 1997                       *****************************/
/********************************************************************************/
#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "own_maths.h"

void unit_vector(double *p_vector);

double size_vector(double *p_vector);

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void rotate_vecs(double *p_vector, double *p_axis, double theta, int num_vecs);

void make_move(atom *p_molecule, atom *p_trial_molecule, int num_host_atoms,
               int num_atoms, double *p_delta, double alpha, int have_host)
{
#include "header.h"

int iloop, index_atom;
double *p_this_delta;
double theta,origin[3], axis[3], vec1[3];
double size;

atom *p_atom, *p_trial, *p_host, *p_trial_host;

/****************************************/
/**** Copy molecule to trial version ****/
/****************************************/


if (alpha > -1E-10 && alpha < 1E-10) return;

if (have_host)
  {
    p_host= p_molecule;
    p_trial_host= p_trial_molecule;

    p_atom= p_molecule+num_host_atoms+1;
    p_trial= p_trial_molecule+num_host_atoms+1;

/*******************************************************************/
/********* Interpret the grad 0 terms as host reverse moves!! ******/
/*******************************************************************/

for (iloop= 0; iloop < num_host_atoms; iloop++)
   {
      p_trial_host->x= p_host->x- alpha* *p_delta;
      p_trial_host->y= p_host->y- alpha* *(p_delta+1);
      p_trial_host->z= p_host->z- alpha* *(p_delta+2);

      p_host++;
      p_trial_host++;
   }
   
/********************************************************************************/
/**** Now align 1->2 with the x-axis ********************************************/
/********************************************************************************/
vec1[0]= *(p_delta+3);
vec1[1]= *(p_delta+4);
vec1[2]= *(p_delta+5);

 if ( vec1[1] != 0 || vec1[2] != 0)
   {
      size= size_vector(&vec1[0]);

      theta= -acos(vec1[0]/size);

      axis[0]=  0;
      axis[1]=  vec1[2]; 
      axis[2]= -vec1[1];

      unit_vector(&axis[0]);

      origin[0]= 0;
      origin[1]= 0;
      origin[2]= 0;

/*DEBUG      rotate(p_molecule, &axis[0], &origin[0], theta, num_host_atoms); */

/********************************************************************************/
/*** If need be do lattice vectors too ******************************************/
/********************************************************************************/

     if (pbc)
       {
/*DEBUG            rotate_vecs(&latt_vec[0], &axis[0], theta, 3); */
        }
  }
/********************************************************************************/
/**** Rotate so that atom three has no z-component ******************************/
/********************************************************************************/
vec1[0]= *(p_delta+6);
vec1[1]= *(p_delta+7);
vec1[2]= *(p_delta+8);

if (vec1[2] != 0)
  {
      size= sqrt(vec1[1]*vec1[1]+vec1[2]*vec1[2]);

      theta= pi-acos(vec1[1]/size);

/************************************************/
/**** above and below y-axis are not       ******/
/**** distuingished use z-co-ord to decide ******/
/**** if a clockwise or anti cw is needed  ******/
/************************************************/

     if (vec1[2] < 0) theta= -theta;
 
     axis[0]= 1; 
     axis[1]= 0;
     axis[2]= 0;

     printf("Atom 3 vector= %10.6f %10.6f %10.6f \n",vec1[0], vec1[1], vec1[2]);
     printf("Final rotation was by %10.6f degrees\n", theta*RAD_TO_DEG);

     origin[0]= 0;
     origin[1]= 0;
     origin[2]= 0;

/*DEBUG      rotate(p_molecule, &axis[0], &origin[0], theta, num_host_atoms); */

/********************************************************************************/
/*** If need be do lattice vectors too ******************************************/
/********************************************************************************/

     if (pbc)
       {
/*DEBUG            rotate_vecs(&latt_vec[0], &axis[0], theta, 3); */
       }
   }

      p_this_delta= p_delta+3;
      p_trial->x= p_atom->x + alpha* *p_this_delta;
      p_trial++;
      p_atom++;
      p_this_delta+= 3;

      p_trial->x= p_atom->x + alpha* *p_this_delta;
      p_this_delta++;
      p_trial->y= p_atom->y + alpha * *p_this_delta;
      p_trial++;
      p_atom++;
      p_this_delta+= 2;

     for (iloop= 3; iloop < num_atoms; iloop++)
       {
         p_trial->x= p_atom->x + alpha * *p_this_delta;
         p_this_delta++;
         p_trial->y= p_atom->y + alpha * *p_this_delta;
         p_this_delta++;
         p_trial->z= p_atom->z + alpha * *p_this_delta;
         p_this_delta++;
         p_atom++;
         p_trial++;
       }   
  }
else
  {
printf("Wheres the host in make_move...\n");
/****************************************/
/**** make Initial move *****************/
/****************************************/
/**** Assume MOPAC format for atoms *****/
/****************************************/

     p_atom= p_molecule+1+num_host_atoms;
     p_trial= p_trial_molecule+1+num_host_atoms;
     p_this_delta= p_delta;

     p_trial->x= p_atom->x + alpha * *p_this_delta;
     p_atom++;
     p_trial++;
     p_this_delta++;

     p_trial->x= p_atom->x + alpha * *p_this_delta;
     p_this_delta++;
     p_trial->y= p_atom->y + alpha * *p_this_delta;
     p_atom++;
     p_trial++;
     p_this_delta++;

     for (iloop= 3; iloop < num_atoms; iloop++)
       {
         p_trial->x= p_atom->x + alpha * *p_this_delta;
         p_this_delta++;
         p_trial->y= p_atom->y + alpha * *p_this_delta;
         p_this_delta++;
         p_trial->z= p_atom->z + alpha * *p_this_delta;
         p_this_delta++;
         p_atom++;
         p_trial++;
       }   
   }

printf("Made %d changes in move\n",p_this_delta-p_delta);
return;
}
