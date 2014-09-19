/**********************************************************************************/
/**** read_input does just that, it obtains the history file name and any *********/
/**** control parameters. *********************************************************/
/**** Started Sept 99 Dave Willock ************************************************/
/**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "maxima.h" 
#include "structures.h" 
#include "reader.h" 
#include "global_values.h"

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

int find_kind(char *token, int level);

double get_double(int *p_error, int skip_lines);

int get_integer( FILE *input_fp, int skip_lines, int *p_error );

int read_input( FILE *input_fp, char *p_title, char *p_history_file,
                char *p_field_file, char *p_rdf_file, char *p_out_file, 
                char *p_pdb_file, char *p_statis_file, int *p_have_statis,
                analysis *p_anal_flags,  
                int *p_have_field, int *p_num_molecules, int *p_num_atoms, 
                int *p_use_type, int *p_num_to_use, int *p_num_anal_flags,
                double *p_cutoff, 
                int *p_have_config,  int *p_have_pdb, 
                atom *p_bisector_atoms,
                double *p_bisect_disp, 
                double *p_min_corr_time, double *p_max_corr_time, double *p_min_msd_time, 
                double *p_max_msd_time, double *p_timestep, int *p_have_rdf, 
                int *p_have_out, int *p_want_raw, int *p_skip_out, 
                twoD_monitors *p_monit_bnds, 
                int *p_num_std_bnds, monitors *p_monitor_set,  monitors *p_monitor_set_angle,
                monitors *p_monitor_set_dihedral, int *p_num_monit, int *p_num_monit_angle,
                int *p_num_monit_dihedral, int *p_num_limit, int *p_num_limit_angle,
                int *p_num_limit_dihedral,
                char *p_xdatcar_file, int *p_have_history, int *p_have_xdatcar,
                char *p_poscar_file, char *p_potcar_file, int *p_have_poscar,
                int *p_have_potcar, windef *p_winref, int *p_need_draw,
                double *p_rgyr_av_start,
                int *p_elli_start, char *p_elli_atype )
 {
   monitors *p_this_moni, *p_this_moni_angle, *p_this_moni_dihedral;
   int iloop, icat, skip, lower_case, leave_case;
   int end_of_input, bail_out, error;
   int *p_this_type;
   int iii,num, index;

   int current_moni, current_limit, current_bin;
   int current_moni_angle, current_bin_angle, current_limit_angle;
   int current_moni_dihedral, current_bin_dihedral, current_limit_dihedral;
   int   token, second_token, third_token;
   char  *tok, *last_tok;

   end_of_input= FALSE;
   skip= TRUE;
   lower_case= TRUE;
   leave_case= FALSE;
   tok = NULL;
  
   *p_have_field = FALSE;
   *p_have_history = FALSE;
   *p_have_config  = FALSE;
   *p_have_pdb = FALSE;
   *p_have_xdatcar = FALSE;
   *p_have_poscar = FALSE;
   *p_have_potcar = FALSE;
   *p_have_rdf = FALSE;
   *p_need_draw = FALSE;
   p_anal_flags--;
   *p_num_monit = -1;
   *p_num_monit_angle = -1;
   *p_num_monit_dihedral = -1;
   *p_num_limit = -1; 
   *p_num_limit_angle = -1; 
   *p_num_limit_dihedral = -1; 
   current_moni=-1;
   current_moni_angle=-1;
   current_moni_dihedral=-1;
   current_limit=-1;
   current_bin=-1;
   current_bin_angle=-1;
   current_bin_dihedral=-1;
   current_limit_angle=-1;
   current_limit_dihedral=-1;

/*** defaults for ellipse calculations ***/
   *p_elli_start=-1;
   strcpy(p_elli_atype, "ALL");

printf("Arrived in read input\n");
   while (!end_of_input)
    {
      tok = tok_get( input_fp, skip, lower_case ); 
      skip=FALSE;

      if ( tok != NULL )
        {
          if ( strcmp(tok , END_OF_INPUT) == 0) end_of_input= TRUE;
          last_tok= tok;

          token = find_kind(tok, PRIME_DIRECTIVE);

          switch (token)
            {
               case TITLE : while ( tok != NULL) 
                               {
                                 tok=tok_get( input_fp, skip, leave_case );
                 
                                 if ( tok != NULL )
                                   {
                                     strcat( p_title, tok);
                                     strcat( p_title, " ");
                                   }
                               } 
                             break;

              case HISTORY_FILE: tok=tok_get( input_fp, skip, leave_case );
                                if (tok != NULL )
                                   {
                                     strcpy( p_history_file, tok);
                                     *p_have_history=TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for history file\n");
                                     exit(0);
                                   }
                                break;

              case CONFIG_FILE: tok=tok_get( input_fp, skip, leave_case );
                                if (tok != NULL )
                                   {
                                     strcpy( p_history_file, tok);
                                     *p_have_config= TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for config file\n");
                                     exit(0);
                                   }
                                break;

              case READ_PDB_FILE:  tok=tok_get( input_fp, skip, leave_case );
                                   if (tok != NULL )
                                     {
                                       strcpy( p_history_file, tok);
                                       *p_have_pdb= TRUE;
                                     }
                                  else
                                     {
                                       printf ("ERROR : No file name given for pdb file\n");
                                       exit(0);
                                     }
                                  break;

              case FIELD_FILE:  tok=tok_get( input_fp, skip, leave_case );
                                if (tok != NULL )
                                   {
                                     strcpy( p_field_file, tok);
                                     *p_have_field = TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for field file\n");
                                     exit(0);
                                   }
                                break;

              case RDF_FILE:  tok=tok_get( input_fp, skip, leave_case );
                              if (tok != NULL )
                               {
                                 strcpy( p_rdf_file, tok);
                                 *p_have_rdf = TRUE;
                               }
                              else
                               {
                                 printf ("ERROR : No file name given for rdf file\n");
                                 exit(0);
                               }
                              break;

              case OUT_FILE:  tok=tok_get( input_fp, skip, leave_case );
                              if (tok != NULL )
                               {
                                 strcpy( p_out_file, tok);
                                 printf("Read out file name as >>%s<<\n", p_out_file);
                                 *p_have_out = TRUE;
                               }
                              else
                               {
                                 printf ("ERROR : No file name given for output file\n");
                                 exit(0);
                               }
                              break;

              case PDB_FILE:  tok=tok_get( input_fp, skip, leave_case );
                              if (tok != NULL )
                               {
                                 strcpy( p_pdb_file, tok);
                                 p_anal_flags->pdb = TRUE;
                               }
                              else
                               {
                                 printf ("ERROR : No file name given for output file\n");
                                 exit(0);
                               }
                              break;

              case STATIS_FILE:  tok=tok_get( input_fp, skip, leave_case );
                              if (tok != NULL )
                               {
                                 strcpy( p_statis_file, tok);
                                 *p_have_statis = TRUE;
                               }
                              else
                               {
                                 printf ("ERROR : No file name given for STATIS file\n");
                                 exit(0);
                               }
                              break;

              case ANALYSE: ++*p_num_anal_flags;
                            p_anal_flags++;
                            if (*p_num_anal_flags == 10)
                              {
                                printf("Maximum number of analyse requests exceeded, limited to 10\n");
                                exit(0);
                              }
                            break;

              case MOI:      if (*p_num_anal_flags < 0)
                               {
                                 printf("moment_of_interia flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->mo_inertia= TRUE;
                               }
                             break;
             
              case ELLI_CALC:  
                             p_anal_flags->elli_calc=TRUE;
                             printf("Found ellipse flag. Starting ellipse calculations\n");
                             tok=tok_get(input_fp, skip, leave_case);
                             second_token=find_kind(tok, SECONDARY_DIRECTIVE);
                             switch(second_token)
                             {
                              case ELLI_START: printf("Found start directive in ellipse\n");
                                             *p_elli_start=get_int(&error, skip);
                                             printf("Start frame for ellipse calculations: %d\n", 
                                                                                 *p_elli_start); 
                                             break;

                              case ELLI_ATYPE: tok=tok_get(input_fp, skip, leave_case);
                                               if(tok != NULL)
                                                {
                                                 strcpy(p_elli_atype, tok);
                                                 printf("Analysing only %s\n",p_elli_atype);
                                                }          
                                               else
                                                { 
                                                printf("Error in reading atom type\n");
                                                exit(0);
                                                }
                                            break;
                              case ELLI_VOL: p_anal_flags->elli_vol=TRUE;
                                             printf("Found volume directive in ellipse. Normalising number of atoms by volume\n");                           
                                            break;
                            }
                             break;

              case LIMITS: printf("Found word LIMITS, setting min/max bondlengths\n"); 
                           skip=FALSE;
                           tok=tok_get( input_fp, skip, leave_case );
                           printf("Second tok = %s\n",tok);
                           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                           printf("That is token %d\n", second_token);
                           switch (second_token)
                                {
                             case BOND2 : (*p_num_limit)++;
                                           current_limit++;

                                         if (current_limit >= MAX_MONIT)
                                           {
                                            printf("ERROR: Too many monitor limit sets requested, current maximum : %d\n", MAX_MONIT);
                                            exit(0);
                                           }

                                          p_this_moni= p_monitor_set+current_limit;

                                          printf("Getting limits for bond lengths\n");
                                          p_this_moni->min_limit= get_double(&error, skip);
                                          printf("Minimum distance: %10.2f \n",p_monitor_set->min_limit);

                                          p_this_moni->max_limit= get_double(&error, skip);
                                          printf("Maximum distance: %10.2f \n",p_monitor_set->max_limit);

                                          break;
 
                             case ANGLE :  (*p_num_limit_angle)++;
                                         current_limit_angle++;
                                         if (current_limit_angle >= MAX_MONIT)
                                           {
                                            printf("ERROR: Too many monitor limit sets requested, current maximum : %d\n", MAX_MONIT);
                                            exit(0);
                                           }
                                          
                                          p_this_moni_angle= p_monitor_set_angle+current_limit_angle;

                                          printf("Getting limits for angles to monitor\n");
                                          p_this_moni_angle->min_limit= get_double(&error, skip);
                                          printf("Minimum angle to bin: %10.2f \n",p_monitor_set_angle->min_limit);
                                          p_this_moni_angle->max_limit= get_double(&error, skip);
                                          printf("Maximum angle to bin: %10.2f \n",p_monitor_set_angle->max_limit);
                                          
                                          printf("Getting limits for angle bond lengths\n");
                                          p_this_moni_angle->max_limit1= get_double(&error, skip);
                                          printf("Maximum distance for atom pair 1-2: %10.2f \n",p_monitor_set_angle->max_limit1);

                                          p_this_moni_angle->max_limit2= get_double(&error, skip);
                                          printf("Maximum distance for atom pair 2-3: %10.2f \n",p_monitor_set_angle->max_limit2);

                                          break;

                             case TORSION : (*p_num_limit_dihedral)++;
                                            current_limit_dihedral++;
                                            if (current_limit_dihedral >= MAX_MONIT)
                                              {
                                               printf("ERROR: Too many monitor limit sets requested, current maximum : %d\n", MAX_MONIT);
                                               exit(0);
                                              }
                                            
                                             p_this_moni_dihedral= p_monitor_set_dihedral+current_limit_dihedral;

                                             printf("Getting limits for dihedral angles to monitor\n");
                                             p_this_moni_dihedral->min_limit= get_double(&error, skip);
                                             printf("Minimum dihedral to bin: %10.2f \n",p_monitor_set_dihedral->min_limit);
                                             p_this_moni_dihedral->max_limit= get_double(&error, skip);
                                             printf("Maximum dihedral to bin: %10.2f \n",p_monitor_set_dihedral->max_limit);
                                          
                                             printf("Getting limits for dihedral angle bond lengths\n");
                                             p_this_moni_dihedral->max_limit1= get_double(&error, skip);
                                             printf("Maximum distance for atom pair 1-2: %10.2f \n", p_monitor_set_dihedral->max_limit1);

                                             p_this_moni_dihedral->max_limit2= get_double(&error, skip);
                                             printf("Maximum distance for atom pair 2-3: %10.2f \n", p_monitor_set_dihedral->max_limit2);

                                             p_this_moni_dihedral->max_limit3= get_double(&error, skip);
                                             printf("Maximum distance for atom pair 3-4: %10.2f \n", p_monitor_set_dihedral->max_limit3);

                                             break;
                                 }
                             break;

              case BINS : printf("Found words BINS, setting number of bins\n"); 
                          tok=tok_get( input_fp, skip, leave_case );
                          second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                          switch (second_token)
                               {
                             
                            case BOND2 : current_bin++;

                                        if (current_bin  >= MAX_MONIT)
                                          {
                                           printf("ERROR: Too many monitor bin sets requested, current maximum : %d\n", MAX_MONIT);
                                           exit(0);
                                          }

                                           p_this_moni= p_monitor_set+current_bin;
                                           p_this_moni->bins= get_integer(input_fp, skip, &error);
                                           printf("Number of bins: %d\n ", p_monitor_set->bins);
                                           break;

                            case ANGLE :   current_bin_angle++;
                                           printf("Setting number of angle bins for instance %d\n", current_bin_angle);
                                           if (current_bin_angle  >= MAX_MONIT)
                                             {
                                              printf("ERROR: Too many monitor angle bin sets requested, current maximum : %d\n", MAX_MONIT);
                                              exit(0);
                                             }
                                           
                                           p_this_moni_angle= p_monitor_set_angle+current_bin_angle;
                                           p_this_moni_angle->bins= get_integer(input_fp, skip, &error);
                                           printf("Number of angle bins: %d\n ", p_monitor_set_angle->bins);
                                           break;

                            case TORSION : /*printf("This one is a torsion, don't know what to do yet!\n");*/
                                           current_bin_dihedral++;
                                           printf("Setting number of dihedral angle bins for instance %d\n", current_bin_angle);
                                           if (current_bin_angle  >= MAX_MONIT)
                                             {
                                              printf("ERROR: Too many monitor dihedral angle bin sets requested, current maximum : %d\n",
                                                     MAX_MONIT);
                                              exit(0);
                                             }

                                           p_this_moni_dihedral= p_monitor_set_dihedral+current_bin_dihedral;
                                           p_this_moni_dihedral->bins= get_integer(input_fp, skip, &error);
                                           printf("Number of diheral angle bins: %d\n ", p_monitor_set_dihedral->bins);
                                           break;
                                }
                            break;

              case MONITOR: printf("Found word moni : \n");              
                            tok=tok_get( input_fp, skip, leave_case );
                            second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                            switch (second_token)
                               {

                             case BOND2 : (*p_num_monit)++;
                                          current_moni++;

                                          if (current_moni >= MAX_MONIT)
                                            {
                                              printf("ERROR: Too many monitors requested, current maximum : %d\n", MAX_MONIT);
                                              exit(0);
                                            }

                                          p_this_moni= p_monitor_set+current_moni;

                                          printf("Getting monitors set %d\n", *p_num_monit);
                                          tok=tok_get( input_fp, FALSE, leave_case );
                                          strcpy( &(p_this_moni->atom1[0]), tok);
                                          printf("First atom: %s\n ",tok);

                                          tok=tok_get( input_fp, FALSE, leave_case );
                                          strcpy( &(p_this_moni->atom2[0]), tok);
                                          printf("Second atom: %s\n ",tok);
                                          printf("\n");
                                          break;
                                        
                              case ANGLE : printf("This one is an angle,\n");
 
                                          (*p_num_monit_angle)++;
                                          current_moni_angle++;

                                          if (current_moni_angle >= MAX_MONIT)
                                           {
                                             printf("ERROR: Too many monitors requested, current maximum : %d\n", MAX_MONIT);
                                             exit(0);
                                           }

                                         p_this_moni_angle= p_monitor_set_angle+current_moni_angle;

                                         printf("Getting monitors set %d\n", *p_num_monit_angle);
                                         tok=tok_get( input_fp, FALSE, leave_case );
                                         strcpy( &(p_this_moni_angle->atom1[0]), tok);
                                         printf("First atom: %s\n ",tok);

                                         tok=tok_get( input_fp, FALSE, leave_case );
                                         strcpy( &(p_this_moni_angle->atom2[0]), tok);
                                         printf("Second atom: %s\n ",tok);
                                         printf("\n");

                                         tok=tok_get( input_fp, FALSE, leave_case );
                                         strcpy( &(p_this_moni_angle->atom3[0]), tok);
                                         printf("Third atom: %s\n ",tok);
                                         printf("\n");
                                    break;

                               case TORSION : /*printf("This one is a torsion, don't know what to do yet!\n");*/
                                             
                                              printf("This one is a dihedral angle,\n");

                                              (*p_num_monit_dihedral)++;
                                              current_moni_dihedral++;

                                              if (current_moni_dihedral >= MAX_MONIT)
                                                {
                                                 printf("ERROR: Too many monitors requested (%d), current maximum : %d\n",
                                                                  current_moni_dihedral, MAX_MONIT);
                                                 exit(0);
                                                }
                                              
                                               p_this_moni_dihedral= p_monitor_set_dihedral+current_moni_dihedral;

                                               printf("Getting monitors set %d\n", *p_num_monit_dihedral);
                                               tok=tok_get( input_fp, FALSE, leave_case );
                                               strcpy( &(p_this_moni_dihedral->atom1[0]), tok);
                                               printf("First atom: %s\n ",tok);

                                               tok=tok_get( input_fp, FALSE, leave_case );
                                               strcpy( &(p_this_moni_dihedral->atom2[0]), tok);
                                               printf("Second atom: %s\n ",tok);
                                               printf("\n");

                                               tok=tok_get( input_fp, FALSE, leave_case );
                                               strcpy( &(p_this_moni_dihedral->atom3[0]), tok);
                                               printf("Third atom: %s\n ",tok);
                                               printf("\n");

                                               tok=tok_get( input_fp, FALSE, leave_case );
                                               strcpy( &(p_this_moni_dihedral->atom4[0]), tok);
                                               printf("Fourth atom: %s\n ",tok);
                                               printf("\n");
                                
                                      break;
                                }
                               break;
 
              case XDATCAR_FILE: printf("Found word XDATCAR\n");
                                tok=tok_get( input_fp, skip, leave_case );
                                if (tok != NULL )
                                   {
                                     strcpy( p_xdatcar_file, tok);
                                     *p_have_xdatcar=TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for XDATCAR file\n");
                                     exit(0);
                                   }
                                break;

              case POSCAR :  printf("Found word POSCAR\n");
                                tok=tok_get( input_fp, skip, leave_case );
                                if (tok != NULL )
                                   {
                                     strcpy( p_poscar_file, tok);
                                     *p_have_poscar=TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for POSCAR file\n");
                                     exit(0);
                                   }

                            break;

              case POTCAR : printf("Found word POTCAR\n");
                            tok=tok_get( input_fp, skip, leave_case );
                                if (tok != NULL )
                                   {
                                     strcpy( p_potcar_file, tok);
                                     *p_have_potcar=TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for POTCAR file\n");
                                     exit(0);
                                   }

                            break;

              case TRAJ : printf("Found word TRAJ, using this XDATCAR file\n"); break;

              case RAW_TRAJ: *p_want_raw= TRUE; break;
 
              case SKIP_OUT: *p_skip_out= get_integer(input_fp, skip, &error); break;

              case ENG_TOT: p_anal_flags->eng_tot= TRUE; break;

              case TEMP_TOT: p_anal_flags->temp_tot= TRUE; break;

              case ENG_CFG: p_anal_flags->eng_cfg= TRUE; break;

              case ENG_VDW: p_anal_flags->eng_vdw= TRUE; break;

              case ENG_COU: p_anal_flags->eng_cou= TRUE; break;

              case ENG_BND: p_anal_flags->eng_bnd= TRUE; break;

              case ENG_ANG: p_anal_flags->eng_ang= TRUE; break;

              case ENG_DIH: p_anal_flags->eng_dih= TRUE; break;

              case ENG_TET: p_anal_flags->eng_tet= TRUE; break;

              case ENG_PV: p_anal_flags->eng_pv= TRUE; break;

              case TEMP_ROT: p_anal_flags->temp_rot= TRUE; break;

              case VIR_CFG: p_anal_flags->vir_cfg= TRUE; break;

              case VIR_VDW: p_anal_flags->vir_vdw= TRUE; break;

              case VIR_COU: p_anal_flags->vir_cou= TRUE; break;

              case VIR_BND: p_anal_flags->vir_bnd= TRUE; break;

              case VIR_ANG: p_anal_flags->vir_ang= TRUE; break;

              case VIR_CON: p_anal_flags->vir_con= TRUE; break;

              case VIR_TET: p_anal_flags->vir_tet= TRUE; break;

              case VOLUME: p_anal_flags->volume= TRUE; break;

              case TEMP_SHL: p_anal_flags->temp_shl= TRUE; break; 

              case ENG_SHL: p_anal_flags->eng_shl= TRUE; break;

              case VIR_SHL: p_anal_flags->vir_shl= TRUE; break;

              case ALPHA: p_anal_flags->alpha= TRUE; break;

              case BETA: p_anal_flags->beta= TRUE; break;

              case GAMMA: p_anal_flags->gamma= TRUE; break;

              case VIR_PMF: p_anal_flags->vir_pmf= TRUE; break;

              case PRESS: p_anal_flags->press= TRUE; break;

              case COORDS: p_anal_flags->coords= TRUE; break;

              case SUM_COMPS: p_anal_flags->sum_comps= TRUE; break;

              case CORRELATION: p_anal_flags->correlation = TRUE;
                                *p_min_corr_time= get_double(&error, skip);
                                *p_max_corr_time= get_double(&error, skip);
                                break;

              case MSD : p_anal_flags->msd = TRUE;
                         *p_min_msd_time= get_double(&error, skip);
                         *p_max_msd_time= get_double(&error, skip);
                         break;

              case LARGE: if (*p_num_anal_flags < 0)
                               {
                                 printf("largest flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->big_med_small[0]= TRUE;
                               }
                             break;

              case MIDDLE: if (*p_num_anal_flags < 0)
                               {
                                 printf("middle flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->big_med_small[1]= TRUE;
                               }
                             break;

              case SMALLEST: if (*p_num_anal_flags < 0)
                               {
                                 printf("smallest flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->big_med_small[2]= TRUE;
                               }
                             break;

              case ANGLES: if (*p_num_anal_flags < 0)
                               {
                                 printf("smallest flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->angles= TRUE;
                               }
                             break;

              case DIPOLES: if (*p_num_anal_flags < 0)
                               {
                                 printf("dipoles flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->dipoles= TRUE;
                               }
                             break;

              case FORSTER: if (*p_num_anal_flags < 0)
                               {
                                 printf("forster flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->forster= TRUE;
                               }
                             break;

              case BONDS:    p_anal_flags->bonds= TRUE;
     
                             *p_num_std_bnds = -1;
                             printf("Getting bonds to monitor\n");
                             tok=tok_get( input_fp, FALSE, leave_case );

/**** Read in range for distances and angles if tok is second directive, ****/
/**** otherwise read as atom labels                                      ****/

                             if (strncmp(tok, "rang", 4) == 0)
                                {
                                  p_monit_bnds->min_limit1 = get_double(&error, skip);
                                  p_monit_bnds->max_limit1 = get_double(&error, skip);
                                  p_monit_bnds->min_limit2 = get_double(&error, skip);
                                  p_monit_bnds->max_limit2 = get_double(&error, skip);
                                }
                             else if (strncmp(tok, "delt", 4) == 0)
                                {
                                  p_monit_bnds->delta1 = get_double(&error, skip);
                                  p_monit_bnds->delta2 = get_double(&error, skip);
                                }
                             else
                                {
                                      printf("tok=%s\n",tok);
                                      if (tok)
                                        {
                                          strcpy( &(p_monit_bnds->atom1[0]), tok);
                                        }
                                       else
                                        {
                                          printf("ERROR: No atom labels given with bond directive\n");
                                          exit(0);
                                        }

                                      tok=tok_get( input_fp, FALSE, leave_case );
                                      printf("tok=%s\n",tok);
                                      if (tok)
                                        {
                                          strcpy( &(p_monit_bnds->atom2[0]), tok);
                                        }
                                       else
                                        {
                                          printf("ERROR: Only one atom label given with bond directive\n");
                                          exit(0);
                                        }
                                }
                             break;

              case NUM_MOLS: *p_num_molecules= get_integer(input_fp, skip, &error); 
                             if (error)
                                {
                                  printf("ERROR: num_molecules directive given without valid number on same line\n");
                                  exit(0);
                                }
                             break;

              case NUM_ATOMS: *p_num_atoms= get_integer(input_fp, skip, &error); 
                               if (error)
                                {
                                  printf("ERROR: atom_number directive given without valid number on same line\n");
                                  exit(0);
                                }
                              break;

/*************************************************************************************/
/*** give the molecule type indicies to be used in the analysis **********************/
/*************************************************************************************/

              case USE_TYPES: p_this_type= p_use_type;
                              for (iloop=0; iloop < MAXMOL; iloop++) 
                                {
                                  *p_this_type = 0;
                                  p_this_type++;
                                }

                              *p_num_to_use = -1;
                              num= get_integer(input_fp, skip, &error);

                              while (!error)
                                {
                                  *(p_use_type+num-1) = 1;
                                  ++*p_num_to_use;
                                  num= get_integer(input_fp, skip, &error);
                                }
                              break;
      
              case CONTACT: if (*p_num_anal_flags < 0)
                               {
                                 printf("contact flag given before analysis directive\n");
                                 exit(0);
                               }
                             else
                               {
                                 p_anal_flags->contact= TRUE;
                                 *p_cutoff= get_double(&error, skip);
                               }
                             break;

              case TIME_STEP: *p_timestep= get_double(&error, skip);
                              break;

              case ADD_ATOM: p_anal_flags->add_atom = TRUE;
                             tok=tok_get( input_fp, skip, leave_case );
                             second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                             switch (second_token)
                               {
                                 case BISECTOR : p_anal_flags->bisector = TRUE;
                                                 for (iloop=0; iloop < 3; iloop++)
                                                   {
                                                     tok=tok_get( input_fp, skip, leave_case );
                                                     strcpy(p_bisector_atoms->label, tok);
                                                     p_bisector_atoms++;
                                                   }
                                                 *p_bisect_disp= get_double(&error, skip);

                                                 break;
                               }
                             break;

              case RADIUS  : printf("Found radius flag need additional input data\n");
                             tok=tok_get( input_fp, skip, leave_case );
                             second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                             switch (second_token)
                               {
                                 case GYRATION : p_anal_flags->rgyr = TRUE;
                                                 printf("Asked for radius of gyration calculation\n");
                                                 break;
                          
                                 case AVERAGE : printf("Found average directive in radius\n");
                                                *p_rgyr_av_start= get_double(&error, skip);
                                                printf("Start time: %10.6f\n", *p_rgyr_av_start);
                                                break;
                               }
                             break;

              case WINDOWS :  p_anal_flags->windows  = TRUE;
                              p_winref->num_atoms = get_integer(input_fp, skip, &error);
                               if (!error)
                                 {
                                   printf("Found windows flag need additional input data for %d atom definitions\n",
                                                      p_winref->num_atoms );

                                   tok=tok_get( input_fp, FALSE, leave_case );
                                   printf("Additional tok = %s\n", tok);

                                   if (tok && strncmp(tok, "per_", 4) == 0)
                                      {
                                         p_winref->num_windows = get_integer(input_fp, skip, &error);
                                      }
                                   else
                                      {
                                         printf("ERROR: windows directive given without expected number of \n");
                                         printf("ERROR: windows per molecule. Use format:\n");
                                         printf("\nwindows 3 per_mol 4\n");
                                         exit(0);
                                      }
                                 }
                               else
                                 {
                                   printf("ERROR: Window directive given without number of atom definitions to expect.\n");
                                   printf("syntax should be like: windows 3 per_mol 4\n");
                                   printf("for a window definition involving 3 atomsi with 4 windows per molecule.\n");
                                   exit(0);
                                 }

                              for (index=0; index<=p_winref->num_atoms-1 ; index++)
                                {
                                  tok=tok_get( input_fp, TRUE, leave_case );
                                  second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                                  switch (second_token)
                                    {
                                      case ATOM : tok=tok_get( input_fp, skip, leave_case );
 
                                                  if (tok)
                                                    {
                                                      strcpy(&(p_winref->atoms[index][0]), tok);
                                                      printf("Read atom: %s\n", p_winref->atoms[index]);
                                                    }                                                    
                                                  else
                                                    {
                                                      printf("ERROR: No atom defined with atom for window\n");
                                                      printf("syntax should be like: atom hc_2 neigh cp_2\n");
                                                      exit(0);
                                                    }

                                                  tok=tok_get( input_fp, skip, leave_case );
                                                  printf("Then got tok >>%s<<\n", tok);

                                                  if (strncmp(tok, "neigh", 5) == 0)
                                                    {
                                                      tok=tok_get( input_fp, skip, leave_case );
                                                      strcpy(&(p_winref->neigh[index][0]), tok);
                                                    }                                                    
                                                  else
                                                    {
                                                      printf("ERROR: Atom defined without neighbour for window\n");
                                                      printf("syntax should be like: atom hc_2 neigh cp_2\n");
                                                      exit(0);
                                                    }

                                                  break;
                                          
                                        default:  printf("ERROR: reading window definition expecting three atom lines:\n");
                                                  printf("syntax should be like: atom hc_2 neigh cp_2\n");
                                                  exit(0);
                                                  break;
                                    }
                                }
                              for (index=0; index<=p_winref->num_atoms-1; index++)
                                {
                                  printf("Read: atom %s with neigh %s\n",p_winref->atoms[index],
                                                                         p_winref->neigh[index]);
                                }
                              break;
                             

                    case DRAW : *p_need_draw = TRUE;
                                 printf("Found drawing flag need additional input data\n");
                                 tok=tok_get( input_fp, skip, leave_case );
                                 second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                                 switch (second_token)
                                   {
                                     case LINE     : p_anal_flags->drawline = TRUE;
                                                     printf("Asked for line drawings.\n");
                                                     break;
                                                 
                                     case ELLIPSE :  p_anal_flags->drawellipse = TRUE;
                                                     printf("Asked for line ellipse.\n");
                                                     break;
                                   }                
                             break;
            }
        }
      skip= TRUE;
    }
 return 0;
 }
