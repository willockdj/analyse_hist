/********************************************************************************/
/**** line search algorithm for MOPAC minimisations *****************************/
/**** Dave Willock April 1997                       *****************************/
/********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

void make_move(atom *p_molecule, atom *p_trial_molecule, int num_atoms,
               double *p_delta, double alpha, int have_host);

double bra_x_ket( double *p_a, double *p_b, int dimension );

void open_file(FILE **p_file, char *p_filename, char *p_status);

void do_mopac(FILE *file_fp, atom *p_molecule, int num_atoms );

void read_mopac_out(FILE *file_fp, atom *p_molecule, int num_atoms,
                    double *p_heat, double *p_grad, int num_grads, 
                    int have_host );

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *p_templ, int num_t_atoms,
                      int *p_need_grad, double *p_grad);

void print_energy(FILE *fp, energy *p_energy, atom *p_molecule,
                                        int num_atoms, int madalung);

void orientate_for_mopac(atom *p_molecule, int *p_num_mol_members, int mopac_one,
                         int total_atoms, double *p_latt_vec, double *p_recip_latt_vec);

void line_search(atom *p_molecule, atom *p_trial_molecule, int num_host_atoms,
                 int num_atoms, int total_atoms,
                 double *p_delta, int num_delta, double *p_grad2, 
                 double *p_alpha, double gd1, double *p_heat2, int have_host,
                 int *p_need_grad )
{
#include "header.h"

int happy, iloop, index_atom;
int num_neg_curve;
int reversed;
int extrapolating, interpolating;

double *p_this_delta;
double gd1_local, gd2;
double abs_gd1, abs_gd2;
double extra, mag_grad;
double alpha_orig, factor;
double sign;

char mopac_filename[FILELEN_MAX];

atom *p_atom;

FILE *mopac_fp;

happy= FALSE;
reversed= FALSE;
extrapolating= FALSE;
interpolating=FALSE;
num_neg_curve=0;
gd1_local= gd1;
abs_gd1= gd1;
if (abs_gd1 < 0) abs_gd1= -abs_gd1;

extra=0.13;
factor=1.0;
alpha_orig= *p_alpha;

while (!happy)
  {

/****************************************/
/**** make Initial move *****************/
/****************************************/

     printf("About to do move: (alpha= %10.6f)\n",*p_alpha);
     for (iloop=0; iloop < num_delta; iloop++) printf("%10.6f\n",*(p_delta+iloop));

     make_move(p_molecule+num_host_atoms, p_trial_molecule+num_host_atoms, 
                                              num_atoms, p_delta, *p_alpha, have_host); 

     if (have_host) orientate_for_mopac(p_trial_molecule, &num_host_atoms, 1, total_atoms,
                                        &latt_vec[0], &recip_latt_vec[0]);


/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
for (iloop=0; iloop < num_atoms; iloop++)
  {
     printf("%s %10.6f  %10.6f  %10.6f \n",
                      (p_trial_molecule+num_host_atoms+iloop)->label, 
                      (p_trial_molecule+num_host_atoms+iloop)->x, 
                      (p_trial_molecule+num_host_atoms+iloop)->y, 
                      (p_trial_molecule+num_host_atoms+iloop)->z);
  }
if (    (p_trial_molecule+num_host_atoms+2)->z >  1e-6 
     || (p_trial_molecule+num_host_atoms+2)->z < -1e-6) exit(0);

/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
/********************************************************************************/
/**** Launch MOPAC **************************************************************/
/********************************************************************************/

     strcpy(mopac_filename, "mopac.dat");

     open_file(&mopac_fp, &mopac_filename[0], "w");

     do_mopac(mopac_fp, p_trial_molecule+num_host_atoms, num_atoms);

/********************************************************************************/
/**** Read MOPAC out file *******************************************************/
/********************************************************************************/

     strcpy(mopac_filename, "mopac.out");

     open_file(&mopac_fp, &mopac_filename[0], "r");

     read_mopac_out(mopac_fp, p_trial_molecule+num_host_atoms, num_atoms, p_heat2, p_grad2, 
                    num_delta, have_host);

     printf("Heat of formation at trial position = %10.6f\n", *p_heat2);

/********************************************************************************/
/**** Now do the host-guest interaction *****************************************/
/********************************************************************************/

     if (have_host)
        {
           calculate_energy(p_trial_molecule, num_host_atoms-1,
                            p_trial_molecule+num_host_atoms, num_atoms-1,
                            p_need_grad, p_grad2);

           print_energy(stdout, &interaction_energy, 
                        p_trial_molecule+num_host_atoms, num_atoms-1, TRUE);
           
        }

     printf("just before gd2 calc alpha = %10.6f\n", *p_alpha);
     printf("GRAD NOW       DELTA NOW\n");

     for (iloop=0; iloop < num_delta; iloop++)
       {
         printf("%10.6f       %10.6f\n",*(p_grad2+iloop), *(p_delta+iloop));
       }
     gd2= *p_alpha * bra_x_ket(p_grad2, p_delta, num_delta);

     printf("GD1_local= %10.6f    GD2 = %10.6f\n", gd1_local, gd2);

     abs_gd2= gd2;
     if (abs_gd2 < 0) abs_gd2= -abs_gd2;

     if ( abs_gd2 < 0.2*abs_gd1 ) 
       {
          happy = TRUE;

          mag_grad=0;
          for (iloop=0; iloop < num_delta; iloop++) mag_grad+= *(p_grad2+iloop)* *(p_grad2+iloop);
          mag_grad= sqrt(mag_grad);

          printf("Thats fine: GD2= %10.6f mag_grad=  %10.6f\n", gd2, mag_grad);
          printf("           grad                    delta\n");
          for (iloop=0; iloop < num_delta; iloop++)
            {
               printf("%d     %10.6f    %10.6f\n", iloop, *(p_grad2+iloop), *(p_delta+iloop));
            }
       }
    else
       {
/********************************************************************************/
/**** Adjust alpha **************************************************************/
/********************************************************************************/

          if (gd1_local < 0 && gd2 < 0)
            {
               if (gd2 > gd1_local)
                 {
                   if (!interpolating) extra= extra/2;
                   interpolating= TRUE;
                   extrapolating= FALSE;
                   factor += extra;
                   *p_alpha= factor * alpha_orig;
                   printf("Lets : Linear extrapolation : alpha now = %10.6f\n",*p_alpha);
                 }
               else if (gd2 < gd1_local)
                 {
                   if (num_neg_curve == 0) extra= extra/2;
                   factor -= extra; 
                   *p_alpha= factor * alpha_orig;
                   num_neg_curve++;
                   if (num_neg_curve > 10)
                     {
                        printf("Negative Curvatute: This has gone on too stopping!\n");
                        exit(0);
                     }
                   printf("Lets : Negative Curvatute: Continuing: alpha = %10.6f\n", *p_alpha);
                 }

            } 
          else if (gd1_local > 0 && gd2 < 0)
            {
               extra=extra/2; 
               factor -= extra;
               *p_alpha= factor * alpha_orig;
               printf("Lets : GD1 positive and GD2 negative: Decreasing alpha: alpha= %10.6f\n",*p_alpha);
            }
          else if (gd1_local < 0 && gd2 > 0)
            {

               factor -= extra;
               *p_alpha= alpha_orig * factor; 
               printf("Lets : Linear Interpolation alpha now %10.6f\n", *p_alpha);
            }
          else if (gd1_local > 0 && gd2 > 0)
            {
               if (!extrapolating) extra= extra/2;
               interpolating= FALSE;
               extrapolating= TRUE; 
               factor -= extra;
               *p_alpha= factor * alpha_orig;
               printf("Lets : Both GD1 and GD2 are positive, Decreasing alpha: alpha= %10.6f\n", *p_alpha);
            }
/***DEBUG          gd1_local=gd2; ******/
          if (*p_alpha < 1E-6 && *p_alpha > -1E-6) *p_alpha= 1E-6; 
      }
  }

return;
}


