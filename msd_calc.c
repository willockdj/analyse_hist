/**********************************************************************/
/*** Mean square displacement calculated over time origins  ***********/
/*** Dave Willock Sept 04                                   ***********/
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"

void msd_calc(vector *p_list, double *p_msd, double *p_times,
              double time_step, int min_msd_ind_diff, 
              int max_msd_ind_diff, int num_frames, int num_mols )
  {
int itime, last_start, istart, start_this_mol_set;
int end_this_mol_set, imol;
int n_sampled;

vector *p_list_start;
vector *p_list_end;

    printf("DEBUG:: In msd_calc have num_frames= %d and num_mols = %d\n", 
                                                             num_frames, num_mols);
    printf("DEBUG:: minimum index difference between frames requested as %d\n", 
                                                                min_msd_ind_diff);
    printf("DEBUG:: maximum index difference between frames requested as %d\n", 
                                                                max_msd_ind_diff);

    if (max_msd_ind_diff > num_frames)
       {
         printf("ERROR : Cannot carry out dot msdelation with times");
         printf(" longer than frames available\n");
         printf("        Fix the input\n");
         exit(0);
       }

/***************************************/
/***** loop over msdelation times *****/
/***************************************/

    for ( itime = min_msd_ind_diff; itime <= max_msd_ind_diff; itime++) 
       {

          *p_times = itime * time_step;
          last_start = num_frames - itime;
          *p_msd=0.0;
          n_sampled=0;
/***************************************/
/***** loop over time origins **********/ 
/***************************************/ 

          for (istart= 0; istart < last_start; istart++)
             {

               start_this_mol_set = (num_mols+1)*istart;
               end_this_mol_set   = (num_mols+1)*(istart+itime);
 
/***************************************/
/***** loop over molecules *************/
/***************************************/
                for ( imol = 0; imol <= num_mols; imol++) 
                   {
                      p_list_start = p_list + start_this_mol_set + imol;
                      p_list_end   = p_list + end_this_mol_set   + imol;

                      *p_msd += (p_list_start->x - p_list_end->x)*
                                (p_list_start->x - p_list_end->x)
                               +(p_list_start->y - p_list_end->y)*
                                (p_list_start->y - p_list_end->y)
                               +(p_list_start->z - p_list_end->z)*
                                (p_list_start->z - p_list_end->z);
                      n_sampled++;
                   }
             }
          printf("Sampled %d time/mol combinations\n", n_sampled);

          *p_msd = *p_msd /(double) n_sampled;
          p_msd++;
          p_times++;
       }
    return;
  }

