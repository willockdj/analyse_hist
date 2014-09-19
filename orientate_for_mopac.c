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

void gen_rot_matrix(double *p_quarternion, double *p_axis, double theta);

void apply_matrix(atom *p_molecule, double *p_matrix, int num_atoms);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void rotate_vecs(double *p_vector, double *p_axis, double theta, int num_vecs);

void orientate_for_mopac(atom *p_molecule, int num_atoms, double *p_mopac_inverter) 
{
#include "header.h"

int iloop;

double size, size_vec1, vec1[3], axis[3], origin[3];
double matrix1[9], matrix2[9], matrix3[9], cofactor[9], determ;
double theta;

atom *p_atom, *p_atom_loop;

/***** DEBUG DEBUG *****/
atom copy[10];

/***** DEBUG DEBUG *****/

p_atom= p_molecule;

printf("DEBUG>> Passed %d atoms to orient...\n", num_atoms);

for (iloop=0; iloop < num_atoms; iloop++)
   {
      copy[iloop]= *(p_molecule+iloop);
   }

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/


/********************************************************************************/
/**** Now align 1->2 with the x-axis ********************************************/
/********************************************************************************/

printf("DEBUG>> In orientate atom 0 is %s at %10.6f %10.6f %10.6f \n",
                              p_atom->label, p_atom->x, p_atom->y, p_atom->z);

p_atom++;

/***** p_atom now refers to atom 2 ******/

vec1[0]= p_atom->x;
vec1[1]= p_atom->y;
vec1[2]= p_atom->z;

size= size_vector(&vec1[0]);

theta= acos(vec1[0]/size);

axis[0]=  0;
axis[1]=  vec1[2];
axis[2]= -vec1[1];

unit_vector(&axis[0]);

origin[0]= 0;
origin[1]= 0;
origin[2]= 0;

gen_rot_matrix(&matrix1[0], &axis[0], theta);

apply_matrix(p_molecule, &matrix1[0], num_atoms);

/****DEBUG DEBUG *****/

printf ("Effect of first rotation on molecule:\n");

p_atom_loop = p_molecule;
for (iloop=0; iloop < num_atoms; iloop++)
  {
    printf("%s  %10.6f %10.6f %10.6f\n", p_atom_loop->label,
                                         p_atom_loop->x,
                                         p_atom_loop->y,
                                         p_atom_loop->z);
    p_atom_loop++;
  }
/****DEBUG DEBUG *****/

/********************************************************************************/
/**** Rotate so that atom three has no z-component ******************************/
/********************************************************************************/

p_atom++;

/***** p_atom now refers to atom 3 ******/

vec1[0]= p_atom->x;
vec1[1]= p_atom->y;
vec1[2]= p_atom->z;

size_vec1= vec1[1]*vec1[1];

size= sqrt(size_vec1+vec1[2]*vec1[2]);

size_vec1= sqrt(size_vec1);

theta= acos(size_vec1/size);

/****DEBUG DEBUG****/
printf("Angle to rotate= %10.6f\n", theta*RAD_TO_DEG);
printf("Size_vec1= %10.6f,  Size = %10.6f\n",size_vec1,size);
printf("size_vec1/size= %10.6f\n", size_vec1/size);
/****DEBUG DEBUG****/

/************************************************/
/**** above and below xy-plane are not     ******/
/**** distuingished use z-co-ord to decide ******/
/**** if a clockwise or anti cw is needed  ******/
/************************************************/

if (vec1[2] < 0) theta= -theta;

axis[0]= 1;
axis[1]= 0;
axis[2]= 0;

origin[0]= 0;
origin[1]= 0;
origin[2]= 0;

gen_rot_matrix(&matrix2[0], &axis[0], theta);

apply_matrix(p_molecule, &matrix2[0], num_atoms);

/************************************************/
/*** Generate compound matrix *******************/
/************************************************/

matrix3[0] = matrix2[0]* matrix1[0]+ matrix2[1]* matrix1[3]+ matrix2[2]* matrix1[6];
matrix3[1] = matrix2[0]* matrix1[1]+ matrix2[1]* matrix1[4]+ matrix2[2]* matrix1[7];
matrix3[2] = matrix2[0]* matrix1[2]+ matrix2[1]* matrix1[5]+ matrix2[2]* matrix1[8];

matrix3[3] = matrix2[3]* matrix1[0]+ matrix2[4]* matrix1[3]+ matrix2[5]* matrix1[6];
matrix3[4] = matrix2[3]* matrix1[1]+ matrix2[4]* matrix1[4]+ matrix2[5]* matrix1[7];
matrix3[5] = matrix2[3]* matrix1[2]+ matrix2[4]* matrix1[5]+ matrix2[5]* matrix1[8];

matrix3[6] = matrix2[6]* matrix1[0]+ matrix2[7]* matrix1[3]+ matrix2[8]* matrix1[6];
matrix3[7] = matrix2[6]* matrix1[1]+ matrix2[7]* matrix1[4]+ matrix2[8]* matrix1[7];
matrix3[8] = matrix2[6]* matrix1[2]+ matrix2[7]* matrix1[5]+ matrix2[8]* matrix1[8];


/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/

/****Apply compound matrix to the copy ***/

apply_matrix(&copy[0], &matrix3[0], num_atoms);

for (iloop=0; iloop < num_atoms; iloop++)
   {
    printf("%s  %10.6f %10.6f %10.6f\n", copy[iloop].label,
                                         copy[iloop].x,
                                         copy[iloop].y,
                                         copy[iloop].z);
   }

printf("First operation:\n");
printf("%10.6f  %10.6f  %10.6f\n",matrix1[0],matrix1[1],matrix1[2]);
printf("%10.6f  %10.6f  %10.6f\n",matrix1[3],matrix1[4],matrix1[5]);
printf("%10.6f  %10.6f  %10.6f\n",matrix1[6],matrix1[7],matrix1[8]);

printf("\nSecond operation:\n");
printf("%10.6f  %10.6f  %10.6f\n",matrix2[0],matrix2[1],matrix2[2]);
printf("%10.6f  %10.6f  %10.6f\n",matrix2[3],matrix2[4],matrix2[5]);
printf("%10.6f  %10.6f  %10.6f\n",matrix2[6],matrix2[7],matrix2[8]);

printf("\nCompound operation:\n");
printf("%10.6f  %10.6f  %10.6f\n",matrix3[0],matrix3[1],matrix3[2]);
printf("%10.6f  %10.6f  %10.6f\n",matrix3[3],matrix3[4],matrix3[5]);
printf("%10.6f  %10.6f  %10.6f\n",matrix3[6],matrix3[7],matrix3[8]);

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/

/************************************************/
/*** Generate The reverse matrix ****************/
/************************************************/

cofactor[0]= matrix3[4]*matrix3[8] - matrix3[5]*matrix3[7];
cofactor[1]= matrix3[5]*matrix3[6] - matrix3[3]*matrix3[8];
cofactor[2]= matrix3[3]*matrix3[7] - matrix3[4]*matrix3[6];

cofactor[3]= matrix3[2]*matrix3[7] - matrix3[1]*matrix3[8];
cofactor[4]= matrix3[0]*matrix3[8] - matrix3[2]*matrix3[6];
cofactor[5]= matrix3[1]*matrix3[6] - matrix3[0]*matrix3[7];

cofactor[6]= matrix3[1]*matrix3[5] - matrix3[2]*matrix3[4];
cofactor[7]= matrix3[2]*matrix3[3] - matrix3[0]*matrix3[5];
cofactor[8]= matrix3[0]*matrix3[4] - matrix3[1]*matrix3[3];

determ= matrix3[0]*cofactor[0] + matrix3[1]*cofactor[1] + matrix3[2]*cofactor[2];

/*** With transpose included ********************/

*p_mopac_inverter= cofactor[0]/determ;
*(p_mopac_inverter+4)= cofactor[4]/determ;
*(p_mopac_inverter+8)= cofactor[8]/determ;
 
*(p_mopac_inverter+1)= cofactor[3]/determ;
*(p_mopac_inverter+2)= cofactor[6]/determ;
*(p_mopac_inverter+5)= cofactor[7]/determ;

*(p_mopac_inverter+3)= cofactor[1]/determ;
*(p_mopac_inverter+6)= cofactor[2]/determ;
*(p_mopac_inverter+7)= cofactor[5]/determ;

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
printf("\nDeterminant: %10.6f\n",determ);
printf("\ncofactors:\n");
printf("%10.6f  %10.6f  %10.6f\n",cofactor[0],cofactor[1],cofactor[2]);
printf("%10.6f  %10.6f  %10.6f\n",cofactor[3],cofactor[4],cofactor[5]);
printf("%10.6f  %10.6f  %10.6f\n",cofactor[6],cofactor[7],cofactor[8]);

printf("check on inverse operator:\n");

matrix1[0] = *p_mopac_inverter* matrix3[0]+ 
             *(p_mopac_inverter+1)* matrix3[3]+ 
             *(p_mopac_inverter+2)* matrix3[6];

matrix1[1] = *p_mopac_inverter* matrix3[1]+ 
             *(p_mopac_inverter+1)* matrix3[4]+ 
             *(p_mopac_inverter+2)* matrix3[7];

matrix1[2] = *p_mopac_inverter* matrix3[2]+ 
             *(p_mopac_inverter+1)* matrix3[5]+ 
             *(p_mopac_inverter+2)* matrix3[8];

matrix1[3] = *(p_mopac_inverter+3)* matrix3[0]+ 
             *(p_mopac_inverter+4)* matrix3[3]+ 
             *(p_mopac_inverter+5)* matrix3[6];

matrix1[4] = *(p_mopac_inverter+3)* matrix3[1]+ 
             *(p_mopac_inverter+4)* matrix3[4]+ 
             *(p_mopac_inverter+5)* matrix3[7];

matrix1[5] = *(p_mopac_inverter+3)* matrix3[2]+ 
             *(p_mopac_inverter+4)* matrix3[5]+ 
             *(p_mopac_inverter+5)* matrix3[8];

matrix1[6] = *(p_mopac_inverter+6)* matrix3[0]+ 
             *(p_mopac_inverter+7)* matrix3[3]+ 
             *(p_mopac_inverter+8)* matrix3[6];

matrix1[7] = *(p_mopac_inverter+6)* matrix3[1]+ 
             *(p_mopac_inverter+7)* matrix3[4]+ 
             *(p_mopac_inverter+8)* matrix3[7];

matrix1[8] = *(p_mopac_inverter+6)* matrix3[2]+ 
             *(p_mopac_inverter+7)* matrix3[5]+ 
             *(p_mopac_inverter+8)* matrix3[8];

printf("%10.6f  %10.6f  %10.6f\n",matrix1[0],matrix1[1],matrix1[2]);
printf("%10.6f  %10.6f  %10.6f\n",matrix1[3],matrix1[4],matrix1[5]);
printf("%10.6f  %10.6f  %10.6f\n",matrix1[6],matrix1[7],matrix1[8]);

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/

return;

}







