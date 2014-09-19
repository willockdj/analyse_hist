/**********************************************************************/
/*** Correlation function for the dot product of a list of vectors ****/
/*** p_list is a list of 3 component vectors in x,y,z order ***********/
/*** p_dots is the list of dot products to be returned      ***********/
/*** min_corr_ind_diff is the minimum correlation time required *******/
/*** max_corr_ind_diff is the maximum correlation time required *******/
/*** num_frames is the number of sets of vectors in the list    *******/
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"

void dot_correlation(vector *p_list, double *p_dots, double *p_times,
                     double time_step, int min_corr_ind_diff, 
                     int max_corr_ind_diff, int num_frames, int num_mols )
  {
int itime, last_start, istart, start_this_mol_set;
int end_this_mol_set, imol;
int n_sampled;

vector *p_list_start;
vector *p_list_end;

    printf("DEBUG:: In dot_correlation have num_frames= %d and num_mols = %d\n", num_frames, num_mols);
    printf("DEBUG:: minimum index difference between frames requested as %d\n", min_corr_ind_diff);
    printf("DEBUG:: maximum index difference between frames requested as %d\n", max_corr_ind_diff);

    if (max_corr_ind_diff > num_frames)
       {
         printf("ERROR : Cannot carry out dot correlation with times longer than frames available\n");
         printf("        Fix the input\n");
         exit(0);
       }

/***************************************/
/***** loop over correlation times *****/
/***************************************/

    for ( itime = min_corr_ind_diff; itime <= max_corr_ind_diff; itime++) 
       {

          *p_times = itime * time_step;
          last_start = num_frames - itime;
          *p_dots=0.0;
          n_sampled=0;
/***************************************/
/***** loop over time origins **********/ 
/***************************************/ 

          for (istart= 0; istart < last_start; istart++)
             {

               start_this_mol_set = (num_mols+1)*istart;
               end_this_mol_set   = (num_mols+1)*(istart+itime);
 
               if (istart == 0 || istart == last_start)
                     printf("In dot_correlation start at %d end %d\n", start_this_mol_set, end_this_mol_set);

/***************************************/
/***** loop over molecules *************/
/***************************************/
                for ( imol = 0; imol <= num_mols; imol++) 
                   {
                      p_list_start = p_list + start_this_mol_set + imol;
                      p_list_end   = p_list + end_this_mol_set   + imol;
                      *p_dots += p_list_start->x * p_list_end->x
                                +p_list_start->y * p_list_end->y;
                                +p_list_start->z * p_list_end->z;
                      n_sampled++;
                   }
                printf("Current correl: %10.6f\n", *p_dots / (double) n_sampled);
             }
          printf("dot corr report: min: %d  current : %d max %d\n", min_corr_ind_diff, itime, max_corr_ind_diff);
          printf("Sampled %d time/mol combinations\n", n_sampled);
          fflush(stdin);

          *p_dots = *p_dots / (double) n_sampled;
          p_dots++;
          p_times++;
       }
    return;
  }

