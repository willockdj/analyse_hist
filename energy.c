/************************************************************/
/* energy.c                                                 */
/* Calculates the interaction between the pore and template */
/* Mode of calculation determined by input file             */
/*                                                          */
/* Parameters:                                              */
/*  struct pore, struct template                            */
/* Returns:                                                 */
/*  energy of interaction, acceptance flag                  */
/*                                                          */
/* Started DWL 27/11/94                                     */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "maxima.h"
#include "ewald.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"

double timer(time_t reference_t);

void min_image( double *x, double *y, double *z);

void ewald_sum(atom *p_molecule, int num_molecule_atoms,
               atom *p_crystal_cell, int num_cell_atoms,
               double *p_kvecs, double *p_kvec2,
               double *p_gvec2, int num_kvecs,
               double *p_cos_sum, double *p_sin_sum,
               int zero_first, int self);

void isol_coul(atom *p_mol_A, int num_mol_A_atoms,
               atom *p_mol_B, int num_mol_B_atoms, int zero_first, int self);

double non_bond_r_eps_sixth(atom *p_pore, int num_p_atoms,
                            atom *p_template, int num_t_atoms, int *p_need_grad,
                            double *p_grad);

double non_bond_a_b_geometric(atom *p_pore, int num_p_atoms,
                              atom *p_template, int num_t_atoms, int *p_need_grad,
                              double *p_grad);

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *p_templ, int num_t_atoms,
                      double *p_kvecs, double *p_kvec2, 
                      double *p_gvec2, int num_kvecs,
                      double *p_cos_sum, double *p_sin_sum,
                      int *p_need_grad, double *p_grad)
{

#include "header.h"
int i,j, iatom, imol;
int zero_first, self_term;
double dx,dy,dz,r,v_sum, contrib;
double pore_self_energy;

atom *p_patom;
atom *p_tatom;
atom *p_atom;
atom *p_image;

time_t time_now;
double start_time;
double end_time;
float  tarray[2];

interaction_energy.steric = FALSE;
interaction_energy.charges = 0;
interaction_energy.non_bonded = 0;

/********************************************************************************/
/****** From just steric information we can calculate the sum of close **********/
/****** contact distances                                              **********/
/********************************************************************************/

if (steric == TRUE) 
  {
    for (i=0; i<= num_t_atoms; i++)  			/* loop over template atoms */
      {
        p_tatom = p_templ+i;       		/* pointer to current template atom */

        for (j=0; j<= num_p_atoms; j++)			 /* loop over pore atoms */
          {
             p_patom = p_pore+j;            /* pointer to current pore atom */

/********************************************************************************/
/********************** calc separation *****************************************/
/********************************************************************************/

             dx = p_tatom->x - p_patom->x;
             dy = p_tatom->y - p_patom->y;
             dz = p_tatom->z - p_patom->z;

/********************************************************************************/
/************ if this is a periodic pore use the nearest image convention *******/
/********************************************************************************/

            if (pbc)
               {
                   min_image( &dx, &dy, &dz);
               }


            r = dx*dx + dy*dy + dz*dz;
            r = sqrt(r);

/********************************************************************************/
/********* calc sum vderWaals radii *********************************************/
/********************************************************************************/

            v_sum = p_patom->vdw + p_tatom->vdw; 
            v_sum = v_sum * vdw_scale;                /* mult by scaling factor */

            if (r < v_sum) interaction_energy.steric = TRUE; 
          } 
			
      }
   }

/*=============================================================================*/
/******** end steric ***********************************************************/
/*=============================================================================*/
/*=============================================================================*/
/******** calculate non bond ***************************************************/
/*=============================================================================*/
DEBUG= FALSE;
if (non_bonded) 
  {
/*******************************************************************************/
/******** Decide which potential type we are using and what combination ********/
/******** rules                                                         ********/
/*******************************************************************************/

    if (DEBUG)  printf("DB>> Trying to get non-bond energy\n");

    if ( strcmp(pot_info.type, R_EPS) == 0 )
      {
         if (strcmp(pot_info.combination, SIXTH_POWER) == 0 )
           {

              p_atom= p_pore-1;
              for (iatom=0; iatom <= num_p_atoms; iatom++)
                 {
                   p_atom++;
                   p_atom->vdw_energy=0;
                 }

              p_atom= p_templ-1;
              for (iatom=0; iatom <= num_t_atoms; iatom++)
                 {
                   p_atom++;
                   p_atom->vdw_energy=0;
                 }

     if (DEBUG)
       {
          printf("DB>> Energy before doing anything = %10.6f\n",interaction_energy.non_bonded );
          printf("DB>> Molecule 1 first atom : %s has %d members\n", p_pore->label, num_p_atoms);
          printf("DB>> Molecule 2 first atom : %s has %d members\n", p_templ->label, num_t_atoms);
       }

              interaction_energy.non_bonded =  
                          non_bond_r_eps_sixth(p_pore, num_p_atoms, 
                                               p_templ, num_t_atoms,
                                               p_need_grad, p_grad);

     if (DEBUG) printf("DB>> Energy after interacting with the pore %10.6f\n",interaction_energy.non_bonded );

/*************************************************************************/
/****** Interactions with symmetry images ********************************/
/*************************************************************************/

            if (symm_set)
              {
                 p_image = p_templ;
                 for ( imol=0; imol <= num_symm_ops; imol++)
                   {
                     p_image += num_t_atoms+1;

                     interaction_energy.non_bonded +=
                                                non_bond_r_eps_sixth(p_image, num_t_atoms,
                                                                     p_templ, num_t_atoms,
                                                                     p_need_grad, p_grad);
                   }
                 if (DEBUG)
                     {
                         printf("DB>> Energy after pore and own images = %10.6f\n",
                                                          interaction_energy.non_bonded );
                         printf("Found non-bond energy = %10.6f\n",interaction_energy.non_bonded);
                     }
              }

           }
         else
           {
              printf ("\nERROR: Attempt to use unknown combination ");
              printf ("rules for this potential type\n when calculating vdw energy.\n");
              printf ("Potential type is currently: %s, quoted combining rules: %s\n",
                                                   pot_info.type, pot_info.combination);
              printf ("Program terminating, goodbye.\n");
              exit(EXIT_FAILURE);
           }
      }
    else if ( strcmp(pot_info.type, A_B) == 0 )
      {    
         if (strcmp(pot_info.combination, GEOMETRIC) == 0 )
           {
              p_atom= p_pore-1;
              for (iatom=0; iatom <= num_p_atoms; iatom++)
                 {
                   p_atom++;
                   p_atom->vdw_energy=0;
                 }

              p_atom= p_templ-1;
              for (iatom=0; iatom <= num_t_atoms; iatom++)
                 {
                   p_atom++;
                   p_atom->vdw_energy=0;
                 }

     if (DEBUG) printf("DB>> Energy before doing anything = %10.6f\n",interaction_energy.non_bonded );

              interaction_energy.non_bonded =  
                          non_bond_a_b_geometric(p_pore, num_p_atoms, 
                                                 p_templ, num_t_atoms,
                                                 p_need_grad, p_grad);

     if (DEBUG) printf("DB>> Energy after interacting with the pore %10.6f\n",interaction_energy.non_bonded );

/*************************************************************************/
/****** Interactions with symmetry images ********************************/
/*************************************************************************/

            if (symm_set)
              {
                 p_image = p_templ;
                 for ( imol=0; imol <= num_symm_ops; imol++)
                   {
                     p_image += num_t_atoms+1;

                     interaction_energy.non_bonded +=
                                                non_bond_a_b_geometric(p_image, num_t_atoms,
                                                                       p_templ, num_t_atoms,
                                                                       p_need_grad, p_grad);
                   }
                 if (DEBUG)
                   {
                         printf("DB>> Energy after pore and own images = %10.6f\n",
                                                          interaction_energy.non_bonded );
                         printf("Found non-bond energy = %10.6f\n",interaction_energy.non_bonded);
                  }
              }

           }
         else
           {
              printf ("\nERROR: Attempt to use unknown combination ");
              printf ("rules for this potential type\n when calculating vdw energy.\n");
              printf ("Potential type is currently: %s, quoted combining rules: %s\n",
                                                   pot_info.type, pot_info.combination);
              printf ("Program terminating, goodbye.\n");
              exit(EXIT_FAILURE);
           }
      }
    else
      {
         printf ("\nERROR: Attempt to use unknown potential ");
         printf ("type when calculating vdw energy.\n");
         printf ("Potential type is currently : %s\n", pot_info.type);
         printf ("Program terminating, goodbye.\n");
         exit(EXIT_FAILURE);
      }
  }
/*=============================================================================*/
/******** end non_bonded *******************************************************/
/*=============================================================================*/
/*=============================================================================*/
/******** calculate coulomb energy *********************************************/
/*=============================================================================*/

if (charges) 
  {
    if (pbc)
      {
/*******************************************************************************/
/***** Template-> framework terms **********************************************/
/*******************************************************************************/

         zero_first=TRUE;
         self_term= FALSE;

         ewald_sum( p_templ, num_t_atoms, p_pore, num_p_atoms, p_kvecs, p_kvec2,
                    p_gvec2, num_kvecs, p_cos_sum, p_sin_sum, zero_first, self_term); 

/*******************************************************************************/
/***** Sum potentials **********************************************************/
/*******************************************************************************/
       
         for ( iatom=0; iatom <= num_t_atoms; iatom++)
            {
               p_atom = p_templ + iatom;
              
               p_atom->electrostatic_pot = coul_prefactor * (p_atom->electrostatic_pot);
               contrib= (p_atom->part_chge)* (p_atom->electrostatic_pot);

               interaction_energy.charges += contrib;
               
            } 
if (DEBUG) printf("DB>>Coulomb interaction energy: %10.6f\n", interaction_energy.charges);
      }
    else
      {
         printf("Will calculate coulomb energy without Ewald sum\n");
      }
  }
DEBUG= FALSE;
return;
}
