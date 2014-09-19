/********************************************************************************/
/**** Re-orientate molecule according to MOPAC rules ****************************/
/**** Molecule 1 (i.e. the second molecule) will be the MOPAC candidate *********/
/********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void unit_vector(double *p_vector);

double size_vector(double *p_vector);

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void rotate_vecs(double *p_vector, double *p_axis, double theta, int num_vecs);

void mopac_trans_mat(atom *p_molecule, int *p_num_mol_members, int mopac_one,
                         int total_atoms) 
{
#include "header.h"

int iloop, index, start_mol;

int *p_this_mol;

double size, vec1[3], axis[3], origin[3];
double a_cross_b[3], b_cross_c[3], c_cross_a[3];
double theta, dot_set[9];
double matrix1[9], matrix2[9];

atom *p_mopac_mol;

/********************************************************************************/
/**** First matrix aligns 1->2 with the x axis **********************************/
/********************************************************************************/

p_mopac_mol++;

vec1[0]= p_mopac_mol->x;
vec1[1]= p_mopac_mol->y;
vec1[2]= p_mopac_mol->z;

size= size_vector(&vec1[0]);

theta= acos(vec1[0]/size);

axis[0]=  0;
axis[1]=  vec1[2];
axis[2]= -vec1[1];

unit_vector(&axis[0]);

origin[0]= 0;
origin[1]= 0;
origin[2]= 0;

gen_rot_matrix(p_molecule, &axis[0], &origin[0], theta, total_atoms);

/********************************************************************************/
/*** If need be do lattice vectors too ******************************************/
/********************************************************************************/

if (pbc)
  {
     rotate_vecs(&latt_vec[0], &axis[0], theta, 3);
  }

/********************************************************************************/
/**** Rotate so that atom three has no z-component ******************************/
/********************************************************************************/

p_mopac_mol++;

vec1[0]= p_mopac_mol->x;
vec1[1]= p_mopac_mol->y;
vec1[2]= p_mopac_mol->z;

size= sqrt(vec1[1]*vec1[1]+vec1[2]*vec1[2]);

theta= acos(vec1[1]/size);

/************************************************/
/**** above and below y-axis are not       ******/
/**** distuingishea use z-co-ord to decide ******/
/**** if a clockwise or anti cw is needed  ******/
/************************************************/

if (vec1[2] > 0) theta= -theta;

axis[0]= 1;
axis[1]= 0;
axis[2]= 0;

printf("Atom 3 vector= %10.6f %10.6f %10.6f \n",vec1[0], vec1[1], vec1[2]);
printf("Final rotation was by %10.6f degrees\n", theta*RAD_TO_DEG);

origin[0]= 0;
origin[1]= 0;
origin[2]= 0;


printf("Using axis: %10.6f %10.6f %10.6f \n", axis[0], axis[1], axis[2]);

rotate(p_molecule, &axis[0], &origin[0], theta, total_atoms);

/********************************************************************************/
/*** If need be do lattice vectors too ******************************************/
/********************************************************************************/

if (pbc)
  {
     rotate_vecs(&latt_vec[0], &axis[0], theta, 3);

/****** generate reciprocal lattice vectors ********************************/

     vec_cross( &latt_vec[0], &latt_vec[3], &a_cross_b[0]);
     vec_cross( &latt_vec[3], &latt_vec[6], &b_cross_c[0]);
     vec_cross( &latt_vec[6], &latt_vec[0], &c_cross_a[0]);

     recip_latt_vec[0] = b_cross_c[0]/cell_volume;
     recip_latt_vec[1]= b_cross_c[1]/cell_volume;
     recip_latt_vec[2]= b_cross_c[2]/cell_volume;

     recip_latt_vec[3]= c_cross_a[0]/cell_volume;
     recip_latt_vec[4]= c_cross_a[1]/cell_volume;
     recip_latt_vec[5]= c_cross_a[2]/cell_volume;

     recip_latt_vec[6]= a_cross_b[0]/cell_volume;
     recip_latt_vec[7]= a_cross_b[1]/cell_volume;
     recip_latt_vec[8]= a_cross_b[2]/cell_volume;
  }

return;

}







