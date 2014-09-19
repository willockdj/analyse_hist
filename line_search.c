/********************************************************************************/
/**** line search algorithm for MOPAC minimisations                    **********/
/**** This is the adapted version that uses bracketting of the minimum **********/
/**** Dave Willock April 1997                                          **********/
/********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

void move_n_calc(atom *p_molecule, atom *p_trial_molecule, atom *p_mopac_molecule,
                 int num_host_atoms, int num_atoms, int total_atoms,
                 double *p_delta, int num_delta, double *p_grad2,
                 double alpha, double *p_energy, double *p_heat2, int have_host,
                 int *p_need_grad );

void line_search(atom *p_molecule, atom *p_trial_molecule,  atom *p_mopac_molecule,
                 int num_host_atoms, int num_atoms, int total_atoms,
                 double *p_delta, int num_delta, double *p_grad2,
                 double *p_alpha, double energy1, double *p_heat2, int have_host,
                 int *p_need_grad )
{
#include "header.h"

int happy, iloop, index_atom;
int igrad, num_iters;

double *p_this_delta;
double gd1_local, gd2;
double abs_gd1, abs_gd2;
double extra, mag_grad;
double alpha_orig, factor;
double sign;
double alpha_a, alpha_b, alpha_c;
double energy_a, energy_b, energy_c;
double step_ratio, a_to_b;
double mag_a_to_b;
double grad_test[50];

char mopac_filename[FILELEN_MAX];

atom *p_atom;

FILE *mopac_fp;

happy= FALSE;
num_iters=0;
step_ratio= 0.61803;
alpha_orig= *p_alpha;

alpha_a= 0;
alpha_b= *p_alpha;
a_to_b = alpha_orig;
alpha_c= step_ratio* a_to_b;

energy_a= energy1;
energy_b= 0;
energy_c= 0;

/*** Don't need grad for line search ***/
*p_need_grad=FALSE;
                         
/****************************************/
/**** make Initial move *****************/
/****************************************/

/***DEBUG  move_n_calc( p_molecule, p_trial_molecule, num_host_atoms,         *****/
/***DEBUG               num_atoms, total_atoms, p_delta, num_delta, p_grad2,         *****/
/***DEBUG               alpha_b, &energy_b, p_heat2, have_host, p_need_grad);         *****/

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
/**** Numerical test of derivatives in a given direction       ***/

 *p_need_grad=TRUE;                                                     

printf ("DEBUG>> numbers: num_delta= %d, num_host_atoms= %d, num_atoms= %d, total_atoms= %d\n",
                  num_delta, num_host_atoms, num_atoms, total_atoms); 

 for (iloop=0; iloop < num_delta; iloop++) grad_test[iloop]= 0;          

 move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,               
              num_atoms, total_atoms, p_delta, num_delta, &grad_test[0],   
              0, &energy_b, p_heat2, have_host, p_need_grad);       

 *p_need_grad=FALSE;                                                    
                        
 for (igrad=0; igrad < num_delta; igrad++)   
   {                     
 for (iloop=0; iloop < num_delta; iloop++) *(p_delta+iloop)= 0;      
                     
 *(p_delta+igrad)=1;  
 alpha_c= -1E-4;   
               
 move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,       
              num_atoms, total_atoms, p_delta, num_delta, p_grad2, 
              alpha_c, &energy_a, p_heat2, have_host, p_need_grad); 
                        
 alpha_c= 1E-4;                 
 move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,            
              num_atoms, total_atoms, p_delta, num_delta, p_grad2,      
              alpha_c, &energy_b, p_heat2, have_host, p_need_grad);      
              
 printf("energy_a = %15.9f energy_b = %15.9f Numerical Gradient =%10.6f, actual = %10.6f\n", 
               energy_a, energy_b, (energy_b-energy_a)/(2*alpha_c), grad_test[igrad]);   
   }    
         
exit(0); 
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/

move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
             num_atoms, total_atoms, p_delta, num_delta, p_grad2,
             alpha_c, &energy_c, p_heat2, have_host, p_need_grad);

while (!happy)
  {

/*     if  (energy_a - energy_c < 1E-4 && energy_a - energy_c > -1E-4  */
/*       && energy_a - energy_b < 1E-4 && energy_a - energy_b > -1E-4) */
/*       {                                                             */
/*          printf("Its all the same to me.....\n");                   */
/*          exit(0);                                                   */
/*       }                                                             */
/*************************************/
/**** energy_b is smallest        ****/
/*************************************/
printf("Summary Information at top of loop: a_to_b= %10.6f\nalpha_a/energy_a, alpha_c/energy_c,alpha_b/energy_b\n",
               a_to_b);

printf("   %10.6f     %10.6f     %10.6f\n", alpha_a, alpha_c, alpha_b);
printf("   %10.6f     %10.6f     %10.6f\n", energy_a, energy_c, energy_b);

     if (energy_b < energy_a && energy_b < energy_c)
        {
/***** Extrapolate from b ************/

           alpha_a  = alpha_c;
           energy_a = energy_c;
           alpha_b  = alpha_a+ a_to_b;
           alpha_c  = alpha_a+ step_ratio* a_to_b;

/*************************************/
/*** re-evaluate for b and c *********/
/*************************************/

           move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
                        num_atoms, total_atoms, p_delta, num_delta, p_grad2,
                        alpha_b, &energy_b, p_heat2, have_host, p_need_grad);

           move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
                        num_atoms, total_atoms, p_delta, num_delta, p_grad2,
                        alpha_c, &energy_c, p_heat2, have_host, p_need_grad);

        }
/*************************************/
/**** energy_a is smallest        ****/
/*************************************/

     else if (energy_a < energy_b && energy_a < energy_c)
       {
/*************************************/
/**** Extrapoplate from a ************/
/*************************************/

           alpha_b  = alpha_c;
           energy_b = energy_c;
           alpha_a  = alpha_b- a_to_b;
           alpha_c  = alpha_a+ step_ratio* a_to_b;

/*************************************/
/**** Eval. a and c ******************/
/*************************************/
 
           move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
                        num_atoms, total_atoms, p_delta, num_delta, p_grad2,
                        alpha_a, &energy_a, p_heat2, have_host, p_need_grad);

           move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
                        num_atoms, total_atoms, p_delta, num_delta, p_grad2,
                        alpha_c, &energy_c, p_heat2, have_host, p_need_grad);
                         
       }
/******************************************************/
/**** energy_c is smallest : correctly bracketted *****/
/******************************************************/

     else if (energy_c <= energy_a && energy_c <= energy_b)
       {
/******************************************************/
/**** move in the largest reset a_to_b ****************/
/******************************************************/

          if (energy_a > energy_b)
           {
              alpha_a= alpha_c;
              energy_a= energy_c;
              a_to_b= alpha_b - alpha_a;
              alpha_c= alpha_a+step_ratio*a_to_b;

           }
          else
           {
              alpha_b= alpha_c;
              energy_b= energy_c;
              a_to_b= alpha_b - alpha_a;
              alpha_c= alpha_a+step_ratio*a_to_b;
           }

/******************************************************/
/**** Re-do c *****************************************/
/******************************************************/

          move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
                       num_atoms, total_atoms, p_delta, num_delta, p_grad2,
                       alpha_c, &energy_c, p_heat2, have_host, p_need_grad);
                         
       }

/*** DEBUG ****/
num_iters++;
mag_a_to_b= a_to_b;
if (mag_a_to_b < 0 ) mag_a_to_b= -mag_a_to_b;

/*******************************************************/
/**** Happy if interval is reduced and we are **********/
/**** properly bracketted                     **********/
/*******************************************************/
happy = happy || ( mag_a_to_b < 0.005*alpha_orig && energy_c <= energy_a && energy_c <= energy_b); 

  }

*p_alpha= alpha_c;

/*******************************************************/
/**** One run to get gradients *************************/
/*******************************************************/
*p_need_grad= TRUE;

move_n_calc( p_molecule, p_trial_molecule, p_mopac_molecule, num_host_atoms,
             num_atoms, total_atoms, p_delta, num_delta, p_grad2,
             alpha_c, &energy_c, p_heat2, have_host, p_need_grad);

return;
}

