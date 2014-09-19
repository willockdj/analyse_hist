/********************************************************************************/
/**** This routine performs a move and calculates the energy / grad    **********/
/**** at the new position                                              **********/
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

void make_move(atom *p_molecule, atom *p_trial_molecule, 
               int num_host_atoms, int num_atoms,
               double *p_delta, double alpha, int have_host);

double bra_x_ket( double *p_a, double *p_b, int dimension );

void open_file(FILE **p_file, char *p_filename, char *p_status);

void do_mopac(FILE *file_fp, atom *p_molecule, int num_atoms, int need_grad );

void orientate_for_mopac(atom *p_molecule, int num_atoms, double *p_mopac_inverter);

void read_mopac_out(FILE *file_fp, atom *p_molecule, int num_atoms,
                    double *p_heat, double *p_grad, double *p_reverse_mopac,
                    int num_grads, int have_host, int need_grads );

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *p_templ, int num_t_atoms,
                      int *p_need_grad, double *p_grad);

void print_energy(FILE *fp, energy *p_energy, atom *p_molecule,
                                        int num_atoms, int madalung);

void move_n_calc(atom *p_molecule, atom *p_trial_molecule, atom *p_mopac_molecule,
                 int num_host_atoms, int num_atoms, int total_atoms,
                 double *p_delta, int num_delta, double *p_grad2, 
                 double alpha, double *p_energy, double *p_heat2, int have_host,
                 int *p_need_grad )
{
#include "header.h"

int happy, iloop, index_atom;
int num_iters;
int icomp, igrad;

double *p_this_delta;
double reverse_mopac[9];

char mopac_filename[FILELEN_MAX];

atom *p_atom;

div_t fract;

FILE *mopac_fp;

/****************************************/
/**** Set trial to molecule for start ***/
/****************************************/
printf("DEBUG>> In move_n_calc need_grad= %d\n", *p_need_grad);

for (iloop=0; iloop < total_atoms; iloop++)
                                 *(p_trial_molecule+iloop)= *(p_molecule+iloop);


/****************************************/
/**** make move *************************/
/****************************************/

printf("About to do move: (alpha= %10.6f)\n",alpha);
printf("DELTA ==>>\n");
for (iloop=0; iloop < num_delta; iloop++) printf("%10.6f\n",*(p_delta+iloop));

make_move(p_molecule, p_trial_molecule, num_host_atoms,
                                         num_atoms, p_delta, alpha, have_host); 

/********************************************************************************/
/**** Do the host-guest interaction first ***************************************/
/********************************************************************************/

if (have_host)
   {
      calculate_energy(p_trial_molecule, num_host_atoms-1,
                       p_trial_molecule+num_host_atoms, num_atoms-1,
                       p_need_grad, p_grad2);

      print_energy(stdout, &interaction_energy, 
                   p_trial_molecule+num_host_atoms, num_atoms-1, TRUE);
           
      *p_energy = interaction_energy.non_bonded;

   }

/********************************************************************************/
/**** Launch MOPAC **************************************************************/
/********************************************************************************/
for (iloop=0; iloop < num_atoms; iloop++) *(p_mopac_molecule+iloop)= *(p_trial_molecule+num_host_atoms+iloop);

for (iloop=0; iloop < num_atoms; iloop++)
  {
    printf("%s  %10.6f %10.6f %10.6f\n",(p_mopac_molecule+iloop)->label,
                                        (p_mopac_molecule+iloop)->x,
                                        (p_mopac_molecule+iloop)->y,
                                        (p_mopac_molecule+iloop)->z);
  }
orientate_for_mopac(p_mopac_molecule, num_atoms, &reverse_mopac[0]);

strcpy(mopac_filename, "mopac.dat");

open_file(&mopac_fp, &mopac_filename[0], "w");

do_mopac(mopac_fp, p_mopac_molecule, num_atoms, *p_need_grad); 

/********************************************************************************/
/**** Read MOPAC out file *******************************************************/
/********************************************************************************/

strcpy(mopac_filename, "mopac.out");

open_file(&mopac_fp, &mopac_filename[0], "r");

read_mopac_out(mopac_fp, p_trial_molecule+num_host_atoms, num_atoms, p_heat2, p_grad2, 
               &reverse_mopac[0], num_delta, have_host, *p_need_grad);   

*p_energy+= *p_heat2;
printf("Heat of formation at trial position = %10.6f\n", *p_heat2);

return;
}

