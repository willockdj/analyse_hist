#include <stdio.h>
#include "maxima.h"
#include "structures.h"

/* protype list for this routine */
void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

void cube_roots(double *p_coeffs, double *p_roots);

/************************/
/**** DEBUG Routines ****/

void centre_of_bscat(double *p_c_of_b, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol );

/**** DEBUG Routines ****/
/************************/
/*---------------------------------------------------------------------------*/

/* work out the moment of interia matrix for the molecule based on SANs bscat parameters */

void moments_of_inertia_bscat(atom *p_molecule, int num_atoms, double *p_c_of_b,
                              double *p_m_of_inertia_b, double *p_eigenvals_bscat )
{
  int icomp,this_atom;
  atom *p_atom;

  double bfact;
  double *p_comp0, *p_comp1, *p_comp2, *p_comp3, *p_comp4, *p_comp5, r2;
  double vec[3];

  double cubic_coeffs[3];

/*************************/
/**** DEBUG Variables ****/

  double c_of_b_dum[3];
  double totb_dum;

/**** DEBUG Variables ****/
/*************************/
/*----------------------------------------------------------------------------*/

/****** Zero moment of inertia matrix *******/

   p_comp0= p_m_of_inertia_b;
   p_comp1= p_m_of_inertia_b+1;
   p_comp2= p_m_of_inertia_b+2;
   p_comp3= p_m_of_inertia_b+3;
   p_comp4= p_m_of_inertia_b+4;
   p_comp5= p_m_of_inertia_b+5;

   *p_comp0 =0.0;
   *p_comp1 =0.0;
   *p_comp2 =0.0;
   *p_comp3 =0.0;
   *p_comp4 =0.0;
   *p_comp5 =0.0;

/****** move molecule to centre of mass *****/
/*  printf("Moving molecule to centre of mass at %10.6f %10.6f %10.6f\n", 
                                                *p_c_of_b, *(p_c_of_b+1), *(p_c_of_b+2)); */
  vec[0]= - *p_c_of_b;
  vec[1]= - *(p_c_of_b+1);
  vec[2]= - *(p_c_of_b+2);
  move_molecule(p_molecule, num_atoms, &vec[0]);

/********************************************/
/****** debug *******************************/

  centre_of_bscat(&c_of_b_dum[0], &totb_dum, p_molecule, num_atoms, -1 );
  printf("Centre of bscat for shifted molecule:  %10.6f %10.6f %10.6f\n",
                           c_of_b_dum[0], c_of_b_dum[1], c_of_b_dum[2]);  
  
/****** debug *******************************/
/********************************************/

/****** Calculate M_of_I matrix *************/

   p_atom= p_molecule;
   for (this_atom=0; this_atom <= num_atoms; this_atom++)
     {
       bfact = p_atom->bscat;
/*       printf("Atomic mass of %s (elem >>%s<< = %10.6f ", p_atom->label, 
                                                                 p_atom->elem, atomic_mass);     */
/*       printf("position: %10.6f %10.6f %10.6f\n", p_atom->x, p_atom->y, p_atom->z);            */

       r2= p_atom->x * p_atom->x + p_atom->y * p_atom->y + p_atom->z * p_atom->z;

/******************************************************************************/
/*** Store as a triangular array since the moment of inertia matrix is ********/
/*** symmetric so you can use the elements :                           ********/
/***                                                                   ********/
/***             0  1  2                                               ********/
/***                3  4                                               ********/
/***                   5                                               ********/
/******************************************************************************/
  
       (*p_comp0) += bfact * (r2 - p_atom->x * p_atom->x);           /* 0 */
  
       (*p_comp1) -= bfact * p_atom->x * p_atom->y;                  /* 1 */
  
       (*p_comp2) -= bfact * p_atom->x * p_atom->z;                  /* 2 */
  
       (*p_comp3) += bfact * (r2 - p_atom->y * p_atom->y);           /* 3 */
  
       (*p_comp4) -= bfact * p_atom->y * p_atom->z;                  /* 4 */
  
       (*p_comp5) += bfact * (r2 - p_atom->z * p_atom->z );          /* 5 */

       p_atom++;
     }

/****** move molecule back from centre of mass *****/

   move_molecule(p_molecule, num_atoms, p_c_of_b); 

/****** Work out eigen vectors of matrix, easy in this case as it is 3x3 so *******/
/****** use standard form from Handbook                                     *******/
/****** Assumes cubic term coefficient is 1                                 *******/
/****** So cubic_coeffs[2]= sq coeff, [1]= linear and [0]= constant term    *******/

/*************************/
/****** DEBUG DEBUG ******/
  
/*  *p_comp0=  1.0; */
/*  *p_comp1=  1.0; */
/*  *p_comp2=  2.0; */
/*  *p_comp3=  1.0; */
/*  *p_comp4= -1.0; */
/*  *p_comp5=  1.0; */
     
/****** DEBUG DEBUG ******/
/*************************/

  printf("DEBUG>> Matrix formed in moments_of_inertia_bscat.c : \n");         
  printf("%10.6f %10.6f %10.6f\n", *p_comp0, *p_comp1, *p_comp2);
  printf("%10.6f %10.6f %10.6f\n", *p_comp1, *p_comp3, *p_comp4);
  printf("%10.6f %10.6f %10.6f\n", *p_comp2, *p_comp4, *p_comp5);

  cubic_coeffs[2] = -*p_comp0 - *p_comp3 - *p_comp5;

  cubic_coeffs[1] = -*p_comp4 * *p_comp4 - *p_comp1 * *p_comp1 
                                         - *p_comp2 * *p_comp2;

  cubic_coeffs[1]+= *p_comp0 * *p_comp3 + *p_comp0 * *p_comp5 
                                        + *p_comp3 * *p_comp5;

  cubic_coeffs[0] =      -*p_comp0 * *p_comp3 * *p_comp5  
                    -2.0* *p_comp1 * *p_comp2 * *p_comp4
                       +  *p_comp0 * *p_comp4 * *p_comp4
                       +  *p_comp3 * *p_comp2 * *p_comp2
                       +  *p_comp5 * *p_comp1 * *p_comp1;

  printf("\nGives cubic coeffs:  %10.6f %10.6f %10.6f\n", 
                  cubic_coeffs[0], cubic_coeffs[1], cubic_coeffs[2]); 
/**** Test roots ***/
/*   cubic_coeffs[0]=-56.0; */
/*   cubic_coeffs[1]=-22.0; */
/*   cubic_coeffs[2]=  5.0; */
/*******************/
  cube_roots(&cubic_coeffs[0], p_eigenvals_bscat);

/*  printf("cubic equation coefficients 0= %10.6f 1=%10.6f 2=%10.6f\n",
                                  cubic_coeffs[0], cubic_coeffs[1], cubic_coeffs[2]); */

 printf("eigenvalues of moment of inertia tensor are : %10.6f %10.6f %10.6f\n",
                                                       *p_eigenvals_bscat, 
                                                       *(p_eigenvals_bscat+1), 
                                                       *(p_eigenvals_bscat+2)); 

return;
}
