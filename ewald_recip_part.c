/****************************************************************************************/
/***** Reciprocal space part of Ewald Sum ***********************************************/
/***** Started 12 June 1996 Dave Willock  ***********************************************/
/****************************************************************************************/


#include <stdio.h>
#include <math.h>
#include "structures.h"
#include "maxima.h"
#include "own_maths.h"
#include "ewald.h"

double ewald_recip_part(double *p_kvecs,double *p_kvec2, double *p_gvec2, 
                        int num_kvecs, atom *p_mol_atom, 
                        double *p_cos_sum, double *p_sin_sum, int self)
{
#include "header.h"

int iatom, mol_index, index;
int num_needed[3];

double recip_energy, g_sqrd, k_sqrd, kx, ky, kz;
double rx, ry, rz, k_dot_r;
double pre_factor, arg;
double g_dot_r_mol, g_dot_r_cell;

atom *p_this_cell_atom;

recip_energy=0;

/*********************************************************************/
/********* Loop over reciprocal space vectors ************************/
/*********************************************************************/

for ( index= 0; index <= num_kvecs; index++)
  {
            
     kx= *p_kvecs;
     p_kvecs++;
     ky= *p_kvecs;
     p_kvecs++;
     kz= *p_kvecs;
     p_kvecs++;

     k_sqrd= *p_kvec2; 
     p_kvec2++;
     g_sqrd= *p_gvec2;
     p_gvec2++;

     g_dot_r_mol = two_pi*( kx* (p_mol_atom->x)
                          + ky* (p_mol_atom->y)
                          + kz* (p_mol_atom->z));

/*********************************************************************/
/**********                                                       ****/
/********** prefactor in recip sum = exp ( -G^2 / (4 * kappa^2) ) ****/
/**********                          ---------------------------- ****/
/**********                                    G^2                ****/
/**********                                                       ****/
/********** We have G=2*pi*k; since we define reciprocal vectors  ****/
/**********                   without the 2*pi for wavevectors    ****/
/********** so in terms of k the exp arguement= -pi^2 * k^2 / kappa^2*/
/**********                                                       ****/
/*********************************************************************/

     pre_factor= 1/g_sqrd;

     arg= -k_sqrd*pi_sqrd_over_kappa_sqrd;
     pre_factor= pre_factor*exp(arg);
          
     recip_energy += pre_factor * (cos(g_dot_r_mol) * (*p_cos_sum)
                                  + sin(g_dot_r_mol) * (*p_sin_sum));
     p_cos_sum++; 
     p_sin_sum++;
  }

if (!self)
  {
    return four_pi_over_vol*recip_energy;
  }
else
  {
    return 0.5*four_pi_over_vol*recip_energy;
  }
}
