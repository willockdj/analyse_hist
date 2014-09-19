/*************************************************/
/* C Program to Print Output from a HISTORY file */
/* Dave Willock 2000                             */
/*************************************************/
/* Latest alterations and capabilities           */
/* December 2000 adding Correlation functions    */
/* starting with moi.                            */
/*************************************************/
/* Update for adsorbates in zeolites fuller      */
/* analysis begun Sept 04, Dave Willock.         */
/*                                               */
/*************************************************/
/*************************************************/
/*************************************************/

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include "maxima.h"
# include "structures.h"
# include "global_values.h"
# include "constants.h"

#define MAIN
# include "own_maths.h"
# include "ewald.h"
# include "header.h"
# include "data.h"
#undef MAIN

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void unit_vector(double *p_vector);

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms );

void moments_of_inertia(atom *p_molecule, int num_atoms, double *p_c_of_m,
                           double *p_m_of_inertia, double *p_eigenvals );

void eigen_vec_3b3( double *p_matrix, double eigen_val, double *p_eigen_vec );

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, 
                double *p_recip_latt_vec, coord_flags *p_fix_flags); 

void write_pdb(FILE *pdb_fp, atom *p_molecule, double *p_abc, int num_atoms);

void msd_calc(vector *p_list, double *p_msd, double *p_times,
              double time_step, int min_corr_ind_diff, 
              int max_corr_ind_diff, int num_frames, int num_mols );

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

void mat_transform( double *p_matrix, double *p_u );

void close_contact(atom *p_molecule1, int num_atoms1,
                   atom *p_molecule2, int num_atoms2,
                   double cutoff2, int pbc, int *p_num_found,
                   pair_list *p_clash_list );

int read_input( FILE *input_fp, char *p_title, char *p_history_file,
                char *p_field_file, char *p_rdf_file, char *p_out_file, 
                char *p_pdb_file, char *p_statis_file, int *p_have_statis,
                analysis *p_anal_flags,  int *p_have_field,
                int *p_num_molecules, int *p_num_atoms, int *p_use_type,
                int *p_num_to_use, int *p_num_anal_flags,
                double *p_cutoff, int *p_have_config, atom *p_bisector_atoms,
                double *p_bisect_disp, double *p_min_corr_time, double *p_max_corr_time,
                double *p_min_msd_time, double *p_max_msd_time,
                double *p_timestep, int *p_have_rdf, int *p_have_out,
                int *p_want_raw, int *p_skip_out, 
                types *p_monit_bnds, int *p_num_std_bnds );

void read_field(FILE *file_fp, int *p_num_molecules, int *p_num_this_mol,
                list_partition *p_demarcation, int *p_num_types );

void read_output(FILE *file_fp, double *p_eng_tot, 
                 double *p_temp_tot,  double *p_eng_cfg, 
                 double *p_eng_vdw,   double *p_eng_cou, 
                 double *p_eng_bnd,   double *p_eng_ang, 
                 double *p_eng_dih,   double *p_eng_tet, 
                 double *p_time,      double *p_eng_pv,   
                 double *p_temp_rot,  double *p_vir_cfg,  
                 double *p_vir_vdw,   double *p_vir_cou,   
                 double *p_vir_bnd,   double *p_vir_ang,   
                 double *p_vir_con,   double *p_vir_tet,   
                 double *p_volume,    double *p_temp_shl,
                 double *p_eng_shl,   double *p_vir_shl,
                 double *p_alpha,     double *p_beta,
                 double *p_gamma,     double *p_vir_pmf,
                 double *p_press,     int *p_num_in_list, 
                 int want_rdf,        rdf *p_rdf,
                 int *p_num_rdfs,     int skip_out,
                 analysis *p_anal_flags );

void read_statis(FILE *file_fp, double *p_enthalpy, 
                 double *p_time, int *p_num );

void write_csv(FILE *fp, char *p_title_x, char *p_title_y, 
               char *p_title_z,
               double *p_x, double *p_y, double *p_z,
               int have_z, int num);

void write_outcsv(double *p_eng_tot, 
                  double *p_temp_tot,  double *p_eng_cfg, 
                  double *p_eng_vdw,   double *p_eng_cou, 
                  double *p_eng_bnd,   double *p_eng_ang, 
                  double *p_eng_dih,   double *p_eng_tet, 
                  double *p_time,      double *p_eng_pv,   
                  double *p_temp_rot,  double *p_vir_cfg,  
                  double *p_vir_vdw,   double *p_vir_cou,   
                  double *p_vir_bnd,   double *p_vir_ang,   
                  double *p_vir_con,   double *p_vir_tet,   
                  double *p_volume,    double *p_temp_shl,
                  double *p_eng_shl,   double *p_vir_shl,
                  double *p_alpha,     double *p_beta,
                  double *p_gamma,     double *p_vir_pmf,
                  double *p_press,     double *p_enthalpy,
                  double *p_stat_time, int num_stat_list,
                  int num_in_list, 
                  int want_rdf,        rdf *p_rdf,
                  int num_rdfs,
                  analysis *p_anal_flags );

int  read_hist_info(FILE *file_fp, char *p_title, int *p_pbc, int *p_levcfg,
                    int imcon, int *p_traj_num );

int  read_hist_frame(FILE *file_fp, atom *p_molecule, int num_atoms, 
                     int *p_nstep, double *p_tstep, int pbc, 
                     double *p_latt_vec, int have_config,
                     int levcfg );

void setup_defaults(void);

void min_image( double *x, double *y, double *z);

int mol_dipole( atom *p_molecule, int num_atoms, vector *p_mol_dipole );

void dot_correlation(vector *p_list, double *p_dots, double *p_times,
                     double time_step, int min_corr_ind_diff, 
                     int max_corr_ind_diff, int num_frames, int num_mols );

void generate_neighbours( atom *p_molecule, int num_atoms, 
                          atom_number *p_types, int *p_num_types,
                          int use_pbc, int first_call);

void gather_molecule(atom *p_molecule, int num_atoms, int which_mol);

int main (int argc, char **argv)
{
FILE *fptr;
FILE *history_fp;
char fname_out[20];
double m_float = 0;
int iloop = 0;
int iatom = 0;
int scan_loop = 0;
char c;
int linecount = 0;
char words[50];
char key[20];
int ikey[20];
int timestep[20];
short num_atoms[20];
float int_timestep[20];
float vectors1[20];
double vectors[20];
char atm_name[20];
float atm_mass[20];
float charge[20];
double velocity[20];
double force[20];
char title_store[50];
int icount, want_raw, done_pdb=FALSE;
int mol_nums[MAX_ATOMS];
short i = 0;
atom molecule[MAX_ATOMS];
atom pdb_frame[MAX_ATOMS];
atom last_frame[MAX_ATOMS];

/*****************************************************/
/* Variables for Mom_inert applied to history files **/
/*****************************************************/

char history_file[80];
char field_file[80];
char out_file[80];
char pdb_file[80];
char car_file[80];
char rdf_file[80];
char statis_file[80];

list_partition demarcation[MAXMOL];

int num_molecules, use_type[MAXMOL], num_to_use;
int num_types, num_this_mol[MAXMOL];
int index, max_atom_index, error;
int start_mol, itype, this_molecule;
int jtype, j_index, j_start_mol, j_this_mol;
int num_cc, abc_form;
int have_config, tot_types;
int have_rdf, have_out;
int num_pdb_atoms, num_pdb_writes;
int super[3], icheck;

FILE *field_fp;
FILE *rdf_fp;
FILE *pdb_fp;
FILE *car_fp;
FILE *out_fp;
FILE *stat_fp;
FILE *msd_fp;
FILE *correlation_fp;

int have_field, levcfg, traj_num;
int nstep, ianal, hist_state; 

double tstep[MAXFRAMES];
double msd[MAXFRAMES];
double time_dots[MAXFRAMES];
double c_of_m[3], total_mass;
double m_of_inertia[6], moi_eigenvals[3], moi_eigenvecs[9];
double angle[3];
double big_eigen, scale, vec[3], vec1[3], temp, abc[6];
double norm[3], perp[3], oh_dist, oh_angl;
double small_eigen, cutoff, cutoff2;
vector mole_dips[MAXMOL];
double dip_size, size;
double kappas[MAXMOL3], tot_kappa;
double dx, dy,dz;

int num_kappas, skip_out;

int big_index, small_index;

pair_list clash_list[MAX_ATOMS];

rdf rdf_list[10];
int num_rdfs;

atom_number atom_types[MAXTYPES];
int num_atom_types;

/**********************************************************************/
/**** map_moi_order will map the order of moi i.e. map_moi_order[0] ***/
/**** contains the biggest, 1 the middle and 3 the smallest index   ***/
/**********************************************************************/

int map_moi_order[3];
int num_of_chars, ichar[LINESIZ], index_new, jloop, good_read;
int header_line[100], title_line[100], date_line[100];
int bis_found, bis_index1, bis_index2, bis_index3;

char c_title_line[100];

atom show_axes[10];
atom bisector_atoms[3];
atom new_atom;
atom temp_atom, temp2_atom;

double bisect_disp;

vector image_disp[MAX_ATOMS];

vector cofm_list[MAX_MSD_LIST];
int num_cofm_list;

analysis anal_flags[10];
int num_anal_flags;

/*****************************************/
/* Variables for totaling moi components */
/*****************************************/

tot_vecs all_mols_moi[3];
double rnum;

/*****************************************/
/* Variables for correlation functions   */
/*****************************************/

double time_step;
double min_msd_time, max_msd_time;
int min_msd_ind_diff, max_msd_ind_diff;
double min_corr_time, max_corr_time;
int min_corr_ind_diff, max_corr_ind_diff;
vector sm_moi_vs_time[MAX_CORR_VECS];
vector me_moi_vs_time[MAX_CORR_VECS];
vector la_moi_vs_time[MAX_CORR_VECS];
vector vel_list[MAX_MSD_LIST];

int frame_index, sm_icorr,me_icorr,la_icorr;
int vel_icorr=0;
double sm_moi_dots[MAXFRAMES];
double me_moi_dots[MAXFRAMES];
double la_moi_dots[MAXFRAMES];
double vel_dots[MAXFRAMES];

double eng_tot[MAX_OUT_LIST];
double enthalpy[MAX_OUT_LIST],   stat_time[MAX_OUT_LIST]; 
double temp_tot[MAX_OUT_LIST],   eng_cfg[MAX_OUT_LIST]; 
double eng_vdw[MAX_OUT_LIST],    eng_cou[MAX_OUT_LIST]; 
double eng_bnd[MAX_OUT_LIST],    eng_ang[MAX_OUT_LIST]; 
double eng_dih[MAX_OUT_LIST],    eng_tet[MAX_OUT_LIST]; 
double time[MAX_OUT_LIST],       eng_pv[MAX_OUT_LIST];   
double temp_rot[MAX_OUT_LIST],   vir_cfg[MAX_OUT_LIST];  
double vir_vdw[MAX_OUT_LIST],    vir_cou[MAX_OUT_LIST];   
double vir_bnd[MAX_OUT_LIST],    vir_ang[MAX_OUT_LIST];   
double vir_con[MAX_OUT_LIST],    vir_tet[MAX_OUT_LIST];   
double volume[MAX_OUT_LIST],     temp_shl[MAX_OUT_LIST];
double eng_shl[MAX_OUT_LIST],    vir_shl[MAX_OUT_LIST];
double alpha[MAX_OUT_LIST],      beta[MAX_OUT_LIST];
double gamma[MAX_OUT_LIST],      vir_pmf[MAX_OUT_LIST];
double press[MAX_OUT_LIST];

int num_in_list, have_statis;
int num_stat_list;

/**********************************/
/*** Variables for bond binning ***/
/**********************************/

double bnd_low,bnd_hgh,bnd_delta, dist;
double guess, centre;

int num_std_bnds;

types monit_bnds[10];

int do_this;
int num_bins, ineigh, this_neigh;
int iii, ibin, iguess, bins[MAX_BINS];

/***********************************************/
/**** Updates for arc file writing, Sept 09 ****/
/***********************************************/
coord_flags fix_flags[MAX_ATOMS]; 
int start_frame;

printf("Starting\n");
/****************************/
/* Some maths constants *****/
/****************************/

one_third= 1.0/3.0;
setup_defaults();

/****************************/
/* Gets the Output Filename */
/****************************/

if (argc != 2)
  {
    printf("ERROR: Wrong number of line arguements given, use syntax: anal_hist input_file\n");
    exit(0);
  }

printf("Input file : %s\n", argv[1]);
strcpy(inputfile, argv[1]);

input_fp = fopen(inputfile, "r");

if (input_fp == NULL )
  {
    printf("Anal_hist ERROR: Error opening file >>%s<<\n", inputfile);
    exit(1);
  }

/***********************************/
/* Set defaults for flags etc ******/
/***********************************/
num_to_use= -1;
num_anal_flags = -1;
have_config= FALSE;
have_statis= FALSE;
have_out= FALSE;
want_raw= FALSE;
skip_out=1;

for (iloop=0; iloop < 10; iloop++)
  {
    anal_flags[iloop].mo_inertia      = FALSE;
    anal_flags[iloop].contact         = FALSE;
    anal_flags[iloop].coords          = FALSE;
    anal_flags[iloop].angles          = FALSE;
    anal_flags[iloop].dipoles         = FALSE;
    anal_flags[iloop].forster         = FALSE;
    anal_flags[iloop].add_atom        = FALSE;
    anal_flags[iloop].bisector        = FALSE;
    anal_flags[iloop].sum_comps       = FALSE;
    anal_flags[iloop].correlation     = FALSE;
    anal_flags[iloop].big_med_small[0]= FALSE;
    anal_flags[iloop].big_med_small[1]= FALSE;
    anal_flags[iloop].big_med_small[2]= FALSE;
    anal_flags[iloop].eng_tot         = FALSE;
    anal_flags[iloop].temp_tot        = FALSE;
    anal_flags[iloop].eng_cfg         = FALSE;
    anal_flags[iloop].eng_vdw         = FALSE;
    anal_flags[iloop].eng_cou         = FALSE;
    anal_flags[iloop].eng_bnd         = FALSE;
    anal_flags[iloop].eng_ang         = FALSE;
    anal_flags[iloop].eng_dih         = FALSE;
    anal_flags[iloop].eng_tet         = FALSE;
    anal_flags[iloop].time            = FALSE;
    anal_flags[iloop].eng_pv          = FALSE;
    anal_flags[iloop].temp_rot        = FALSE;
    anal_flags[iloop].vir_cfg         = FALSE;
    anal_flags[iloop].vir_vdw         = FALSE;
    anal_flags[iloop].vir_cou         = FALSE;
    anal_flags[iloop].vir_bnd         = FALSE;
    anal_flags[iloop].vir_ang         = FALSE;
    anal_flags[iloop].vir_con         = FALSE;
    anal_flags[iloop].vir_tet         = FALSE;
    anal_flags[iloop].volume          = FALSE;
    anal_flags[iloop].temp_shl        = FALSE;
    anal_flags[iloop].eng_shl         = FALSE;
    anal_flags[iloop].vir_shl         = FALSE;
    anal_flags[iloop].alpha           = FALSE;
    anal_flags[iloop].beta            = FALSE;
    anal_flags[iloop].gamma           = FALSE;
    anal_flags[iloop].vir_pmf         = FALSE;
    anal_flags[iloop].press           = FALSE;
    anal_flags[iloop].hin             = FALSE;
    anal_flags[iloop].pdb             = FALSE;
    anal_flags[iloop].msd             = FALSE;
    anal_flags[iloop].bonds           = FALSE;
  }

time_step=0.0;
/***********************************/
/* Read directives from input file */
/***********************************/

good_read= read_input( input_fp, &title[0], &history_file[0],
                       &field_file[0], &rdf_file[0], &out_file[0],
                       &pdb_file[0], &statis_file[0], &have_statis,
                       &anal_flags[0], &have_field,
                       &num_molecules, &num_this_mol[0], &use_type[0],
                       &num_to_use, &num_anal_flags, &cutoff, &have_config,
                       &bisector_atoms[0], &bisect_disp, &min_corr_time,
                       &max_corr_time,
                       &min_msd_time, &max_msd_time, &time_step, 
                       &have_rdf, &have_out, &want_raw, &skip_out,
                       &monit_bnds[0], &num_std_bnds );

printf("Back from read_input with out_file: >>%s<<\n", out_file);

/******************************************/
/*** Read in FIELD file if we have one ****/
/******************************************/

if (have_field)
  {
    if ( (field_fp = fopen (&field_file[0],"r") ) == NULL)
      {
        printf ("ERROR Opening FIELD FILE : %s\n", &field_file[0]);
        return 1;
      }

     printf("Getting data from field file >>%s<<\n", field_file);
     read_field( field_fp, &num_molecules, &num_this_mol[0], &demarcation[0],
                 &num_types);
  }

if (have_rdf)
  {
    if ( (rdf_fp = fopen (&rdf_file[0],"r") ) == NULL)
      {
        printf ("ERROR Opening RDF FILE : %s\n", &rdf_file[0]);
        return 1;
      }

    printf("Converting rdf file %s to seperate csv format files\n");

    

    fclose(rdf_fp);
  }

if (have_out)
  {
    if ( (out_fp = fopen (&out_file[0],"r") ) == NULL)
      {
        printf ("ERROR Opening OUT FILE : %s\n", &out_file[0]);
        return 1;
      }

    printf("Analysing dlpoly OUTPUT file : %s\n", &out_file[0]);

    read_output(out_fp, &eng_tot[0], 
                &temp_tot[0],   &eng_cfg[0], 
                &eng_vdw[0],    &eng_cou[0], 
                &eng_bnd[0],    &eng_ang[0], 
                &eng_dih[0],    &eng_tet[0], 
                &time[0],       &eng_pv[0],   
                &temp_rot[0],   &vir_cfg[0],  
                &vir_vdw[0],    &vir_cou[0],   
                &vir_bnd[0],    &vir_ang[0],   
                &vir_con[0],    &vir_tet[0],   
                &volume[0],     &temp_shl[0],
                &eng_shl[0],    &vir_shl[0],
                &alpha[0],      &beta[0],
                &gamma[0],      &vir_pmf[0],
                &press[0],      &num_in_list, 
                have_rdf,       &rdf_list[0],
                &num_rdfs,      skip_out,
                &anal_flags[0] );

    if (have_statis)
     {
       if ( (stat_fp = fopen (&statis_file[0],"r") ) == NULL)
         {
           printf ("ERROR Opening STATIS file : %s\n", &statis_file[0]);
           return 1;
         }
       printf("Going to read_statis\n");
       read_statis(stat_fp, &enthalpy[0], &stat_time[0], &num_stat_list);
       printf("Back from read_statis with %d list\n", num_stat_list);
     }

    printf("Passing lists to write_outcsv\n");

    write_outcsv( &eng_tot[0], 
                  &temp_tot[0],   &eng_cfg[0], 
                  &eng_vdw[0],    &eng_cou[0], 
                  &eng_bnd[0],    &eng_ang[0], 
                  &eng_dih[0],    &eng_tet[0], 
                  &time[0],       &eng_pv[0],   
                  &temp_rot[0],   &vir_cfg[0],  
                  &vir_vdw[0],    &vir_cou[0],   
                  &vir_bnd[0],    &vir_ang[0],   
                  &vir_con[0],    &vir_tet[0],   
                  &volume[0],     &temp_shl[0],
                  &eng_shl[0],    &vir_shl[0],
                  &alpha[0],      &beta[0],
                  &gamma[0],      &vir_pmf[0],
                  &press[0],      &enthalpy[0],
                  &stat_time[0],  num_stat_list,
                  num_in_list, 
                  have_rdf,       &rdf_list[0],
                  num_rdfs,
                  &anal_flags[0] );

    fclose(out_fp);
  }

/******************************************/
/*** Check that use_type is set or       **/
/*** default to using all molecule types **/
/******************************************/

if (num_to_use < 0)
  {
    num_to_use= num_types;
    for (iloop=0; iloop <= num_types; iloop++) use_type[iloop] = 1;
  }
else if (num_to_use > num_types)
  {
    printf("Error: More types to analyse than have been defined\n");
    exit(0);
  }

/***********************************/
/*** Report back about input *******/
/***********************************/

printf("Title                            : %s\n", &title[0]);

if (!have_config)
  {
    printf("HISTORY file                     : %s\n", &history_file[0]);
  }
else
  {
    printf("CONFIG file                      : %s\n", &history_file[0]);
  }

printf("DLPOLY run used a time between frames of %10.6f ps\n", time_step );

printf("%d Analysis passes requested\n", num_anal_flags+1);

for (iloop=0; iloop <= num_anal_flags; iloop++)
  {
    printf("\nAnalysis pass %d is for:\n", iloop+1);
    if (anal_flags[iloop].mo_inertia) printf("Moment of inertia calculation\n");
    if (anal_flags[iloop].big_med_small[0])
      {
         printf("Largest eigenvalue will be reported\n");
         if (anal_flags[iloop].sum_comps)
              printf("Totals and averages of the unit vector components will also be generated\n");
      }
    if (anal_flags[iloop].big_med_small[1])
      {
         printf("Middle eigenvalue will be reported\n");
         if (anal_flags[iloop].sum_comps)
              printf("Totals and averages of the unit vector components will also be generated\n");
      }
    if (anal_flags[iloop].big_med_small[2])
      {
         printf("Smallest eigenvalue will be reported\n");
         if (anal_flags[iloop].sum_comps)
              printf("Totals and averages of the unit vector components will also be generated\n");
      }
    if (anal_flags[iloop].correlation)
                {

                   if (time_step < 0.0) 
                     {
                        printf("ERROR: No time_frame set in input file, cannot carry out correlation analysis\n");
                        exit(0);
                     }
                   else
                     {
                         printf("Will report correlation functions for moi\n");
                         printf("Min time for correlation %10.6f, max time %10.6f\n",
                                                                      min_corr_time, max_corr_time);
                         min_corr_ind_diff= min_corr_time/time_step; 
                         max_corr_ind_diff= max_corr_time/time_step;

                         printf("Corresponding frame separations : min : %d    max : %d\n", 
                                                                         min_corr_ind_diff, max_corr_ind_diff);
                     }
                }
    if (anal_flags[iloop].msd)
                {

                   if (time_step < 0.0) 
                     {
                        printf("ERROR: No time_frame set in input file, ");
                        printf("cannot carry out msd analysis\n");
                        exit(0);
                     }
                   else
                     {
                         printf("Will report msd function\n");
                         printf("Min time for msd %10.6f, max time %10.6f\n",
                                                          min_msd_time, max_msd_time);
                         min_msd_ind_diff= min_msd_time/time_step; 
                         max_msd_ind_diff= max_msd_time/time_step;
                         printf("Corresponding indicies : %d and %d\n",
                                        min_msd_ind_diff, max_msd_ind_diff);

                         if ( min_msd_time < time_step )
                            {
                               printf("ERROR : Minimum time for msd must be greater than\n");
                               printf("        or equal to the time_frame value or there\n");
                               printf("        are no frames to compare!\n");
                               exit(0);
                            }

                     }
                }
    if (anal_flags[iloop].angles)
                   printf("Will report angles between MOI vector and axes\n");
    if (anal_flags[iloop].dipoles)
                   printf("Will report molecular dipole moments\n");
    if (anal_flags[iloop].forster)
                   printf("Will report forster dipole orientation factor averaged over each frame\n");
    if (anal_flags[iloop].coords)
       {
          printf("Will print molecular co-ordinates for each molecule so that all atoms are\n");
          printf("in a common unit cell.\n");
       }
    if (anal_flags[iloop].contact)
      {
        printf("Will search for all inter-molecular close contacts below %10.6f Angstroms\n",
                cutoff);
        cutoff2= cutoff*cutoff;
      }
    if (anal_flags[iloop].add_atom)
      {
        printf("Will add atoms at the geometric positions defined in the input file\n");
        if (anal_flags[iloop].bisector )
          {
             printf("A new site will be added along the bisector of %s %s %s triplets of atoms", 
                           bisector_atoms[0].label, bisector_atoms[1].label, bisector_atoms[2].label);
             printf(" %10.6f from the central atom into the acute angle\n", bisect_disp);
          }
      }
    if (anal_flags[iloop].pdb)
      {
        printf("Will output trajectory file in pdb format as %s\n", pdb_file);
        if ( want_raw )
          {
             printf("Will print raw pdb data, i.e. with periodic boundaries shown\n");
          }
      }
   
   if (anal_flags[iloop].bonds)
     {
       printf("Will create histogram of bond lengths\n");

       printf("For %d bonds involving only: ");
       for ( iii = 0; iii <= num_std_bnds; iii++) 
                           printf("%s ", monit_bnds[iii].name);
       printf("\n");
/*** define range of bonds to look for ***/

       bnd_low   = 0.5;
       bnd_hgh   = 2.0;
       bnd_delta = 0.01;

       num_bins = (bnd_hgh - bnd_low)/bnd_delta;
       printf("Will require %d bins\n", num_bins); 
       if (num_bins > MAX_BINS)
         {
           printf("ERROR: Current maximum number of bins = %d\n",
                      MAX_BINS);
           exit(0);
         }

      for (ibin = 0; ibin < num_bins; ibin++) bins[ibin]=0;
     }        
  }

if (have_field)
     printf("FIELD   file                     : %s\n", &field_file[0]);

printf("Number of molecules per frame    : %d\n", num_molecules);
index = 0;
for (iloop=0; iloop < num_molecules; iloop++)
  {
     printf("Number of occurances of molecule %d : %d\n", iloop+1, 
                                                     num_this_mol[iloop]);

     for (jloop=0; jloop < num_this_mol[iloop]; jloop++)
       {
         printf("Number of atoms in molecule %d occurance %d : %d start : %d end : %d\n", 
                         iloop+1, jloop+1, demarcation[index].num,
                         demarcation[index].start, demarcation[index].end ); 
         index++;
       }
  }
max_atom_index= demarcation[--index].end;


	printf("\n\nAnalysis will look at %d molecule types: ", num_to_use+1);
	for (iloop=0; iloop < num_molecules; iloop++)
	  {
	    if (use_type[iloop]) printf("%d ", iloop+1);
	  }
	printf("\n");

	printf("The highest expected atom index in a frame is %d\n", max_atom_index);

	/***************************************************/
	/*** Check we have everything and all is well ******/
	/***************************************************/

	error= FALSE;
	if (num_molecules <= 0)
	  {
	    printf("Error: The number of molecules simulated has not been given");
	    printf(" in the input or in a FIELD file\n");
	    error= TRUE;
	  }

	index =0;
	for (iloop=0; iloop < num_molecules; iloop++)
	  {
	     if ( num_this_mol[iloop] <= 0 ) 
	       {
		 printf("Error: Molecule %d has no number of occurances set\n", iloop+1);
		 error= TRUE;
	       }

	     for (jloop=0; jloop < num_this_mol[iloop]; jloop++)
	       {
		 if (demarcation[index].num <= 0)
		   {
		     printf("Error: Molecule %d occurance %d has no atom number set.\n",
									   iloop+1, jloop+1);
		     error= TRUE;
		   }
		 index++;
	       }
	  }
	/***********************************/
	/** DEBUG DEBUG DEBUG DEBUG DEBUG **/
	/*exit(0); */
	/** DEBUG DEBUG DEBUG DEBUG DEBUG **/
	/***********************************/

	if (error) return 1;

	/**********************************/
	/* Opens HISTORY File for Reading */
	/**********************************/

	if ( (history_fp = fopen (&history_file[0],"r") ) == NULL)
	    {
	    printf ("ERROR Opening HISTORY FILE\n");
	    perror ("open");
	    return 1;
	    }

	/*******************************/
	/* Scans Information From File */
	/*******************************/

	hist_state = read_hist_info( history_fp, &title_store[0], 
                                     &pbc, &levcfg, have_config, &traj_num);

/***** 1 signifies HISTORY begins with timestep line ***/
        if ( hist_state == 1)
          {
            fclose(history_fp);
	    history_fp = fopen (&history_file[0],"r");
          }

	if (have_config)
	  {
	    printf (" CONFIG  file title : %s\n",title_store); 
	  }
	else
	  {
	    printf (" HISTORY file title : %s\n",title_store); 
	  }

	printf("Level of atomic information (levcfg)           : %d\n", levcfg);
	printf("Periodic boundary status                       : %d\n", pbc);

	if (!have_config)
	     printf("Number of atoms per frame according to HISTORY : %d\n", traj_num);

	/*********************************************************************************/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/*********************************************************************************/

	sm_icorr = 0;
	me_icorr = 0;
	la_icorr = 0;
	abc_form= FALSE;
	frame_index = -1;
        num_cofm_list= -1;
	while ( read_hist_frame( history_fp, &molecule[0], max_atom_index, &nstep, 
                                 &tstep[frame_index+1], 
			         pbc, &latt_vec[0], have_config, levcfg ) > 0) 
	  {
	      frame_index++;

/** Start of H-bond processing additions ***/
/*put all labels and corresponding distances into arrays              */
/*calculate and print average distances for h-bonding for this frame  */
/*move to next frame and repeat calculation                           */
/*get averages across all frames and calculate global average         */

              double cutoff_distance = 4.0;
              double min_distance = 1.0;
              int num_ho_oh = 0;
              int num_ho_oc = 0;
              int num_ho_f  = 0;                
              int iother;
              int is_oh, is_ho, is_oc, is_f, do_other;
              double x1, x2, y1, y2, z1, z2 = 0;
              double distance = 0;
              printf("DEBUG>> Looking for pairs of atoms in range %10.6f to %10.6f\n",
                                        min_distance, cutoff_distance);                                             

              for (icheck=0; icheck<=max_atom_index; icheck++)                  
                 {
                    is_oh = FALSE;
                    is_ho = FALSE;
                    is_oc = FALSE;
                    is_f  = FALSE;
                    do_other= FALSE;

                if (strcmp(molecule[icheck].label,"oh") == 0)
                  {
                    is_oh = TRUE;
                    do_other= TRUE;
                  }
                else if (strcmp(molecule[icheck].label,"ho") == 0)
                       {
                        is_ho = TRUE;
                        do_other = TRUE;
                       }
                else if (strcmp(molecule[icheck].label,"oc") == 0)
                       {
                        is_oc = TRUE;
                        do_other = TRUE;
                       }
                else if (strcmp(molecule[icheck].label,"f") == 0)
                       {
                        is_f = TRUE;
                        do_other = TRUE;
                       }

                if (do_other)
                  {

                     printf("Checking label %d >>%s<< elem >>%s<< x: %10.6f y: %10.6f z: %10.6f\n", icheck,               
                                                          molecule[icheck].label,
                                                          molecule[icheck].elem,
                                                          molecule[icheck].x,
                                                          molecule[icheck].y,
                                                          molecule[icheck].z );
             /*begin hydrogen bonding calculation*/             

                    x1=  molecule[icheck].x;
                    y1=  molecule[icheck].y;
                    z1=  molecule[icheck].z;
             
                    for (iother=icheck+1; iother<=max_atom_index; iother++)                  
                      {
                        x2=  molecule[iother].x;
                        y2=  molecule[iother].y;
                        z2=  molecule[iother].z;

             /*start with first frame*/ 
             /* for each pair of atoms with appropriate labels */
                
                if ((is_oh && strcmp(molecule[iother].label,"ho") == 0)
                    || (strcmp(molecule[iother].label,"oh") == 0 && strcmp(molecule[icheck].label,"ho") == 0))
                  {
                    distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
                    printf("Found case of oh and ho distance: %10.6f\n",distance);
                     if (min_distance <= distance &&  distance < cutoff_distance)            
                       {
                        printf("%s %s %10.6f\n", molecule[icheck].elem, molecule[iother].elem, distance);
                        num_ho_oh++;
                       }
                  }
                else if ((is_ho && strcmp(molecule[iother].label,"oc") == 0)
                        || (strcmp(molecule[iother].label,"ho") == 0 && strcmp(molecule[icheck].label,"oc") == 0))
                  {
                     distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
                     printf("Found case of oh and oc distance: %10.6f\n", distance);
                      if (min_distance <= distance &&  distance < cutoff_distance)            
                        {
                         printf("%s %s %10.6f\n", molecule[icheck].elem, molecule[iother].elem, distance);
                         num_ho_oc++;
                        }
                  }
                else if ((is_ho && strcmp(molecule[iother].label,"f") == 0)
                        || (strcmp(molecule[iother].label,"ho") == 0 && strcmp(molecule[icheck].label,"f") == 0))
                  {
                      distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
                      printf("Found case of oh and f distance: %10.6f\n", distance);
                       if (min_distance <= distance &&  distance < cutoff_distance)            
                         {
                          printf("%s %s %10.6f\n", molecule[icheck].elem, molecule[iother].elem, distance);
                          num_ho_f++;
                         }
                     }
                   }
                 }
               }

              exit(0);
              /*end hydrogen bonding calculation*/

              if (frame_index >= MAXFRAMES)
                {
                   printf("ERROR : Too many frames in HISTORY file for current version.\n");
                   printf("        Recompile with MAXFRAMES increased.\n");
                   exit(0);
                }
              printf("Processing frame %d\n",frame_index);
	      if (!have_config) printf("frame at timestep %d, %10.6f\n\n", 
                                           nstep, tstep[frame_index]);
	      
	      if (pbc > 0)
		{
		  printf("This is a periodic simulation current lattice vectors are:\n");
		  printf("%10.6f %10.6f %10.6f\n", latt_vec[0], latt_vec[1], latt_vec[2]);
		  printf("%10.6f %10.6f %10.6f\n", latt_vec[3], latt_vec[4], latt_vec[5]);
		  printf("%10.6f %10.6f %10.6f\n", latt_vec[6], latt_vec[7], latt_vec[8]);

		  cart_latt_vecs(&abc[0], &latt_vec[0], &real_latt_sizes[0],
				 &recip_latt_vec[0], &recip_latt_sizes[0],
				 &cell_volume, max_atom_index, abc_form);

		}

              generate_neighbours( &molecule[0],  max_atom_index, 
                                   &atom_types[0], &num_atom_types,
                                   TRUE, TRUE);
 
	      for (ianal=0; ianal <= num_anal_flags; ianal++)
		{
	/*********************************************/
	/*** Set frame averages and counters *********/
	/*********************************************/

		  if (anal_flags[ianal].mo_inertia && anal_flags[ianal].sum_comps)
		    {
		      for (iloop=0; iloop < 3; iloop++)
			{
			  all_mols_moi[iloop].totx = 0.0;
			  all_mols_moi[iloop].toty = 0.0;
			  all_mols_moi[iloop].totz = 0.0;
			  all_mols_moi[iloop].num  = 0;
			}
		    }

		  if (anal_flags[ianal].forster )
		    {
		      num_kappas=0;
		      tot_kappa= 0.0;
		    }

		  printf("Carrying out Analysis pass number %d for %d molecules\n", ianal+1, num_molecules);

		  index =0;
		  tot_types=0;
		  for (itype=0; itype < num_molecules; itype++)
		    {
		      for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
			{
			   tot_types++;
			}
		    }

		  for (itype=0; itype < num_molecules; itype++)
		    {
	/*********************************************/
	/*** Only include flagged molecule types   ***/
	/*********************************************/
		      if (use_type[itype])
			 {
                           printf("Analysing molecule type %d\n",itype);

	/*********************************************/
	/*** Consider each of this molecule type *****/
	/*********************************************/

			   for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
			      {
				start_mol= demarcation[index].start;

	/***************************************************************/
	/*** Move all atoms to same unit cell as the first in the list */
	/*** Remember all atom displacements for fixing later          */
	/***************************************************************/
/*
                                if (pbc > 0)
				   {
				     for (iloop=1; iloop <= demarcation[index].num; iloop++)
					{
					  image_disp[iloop].x= molecule[start_mol].x- molecule[start_mol+iloop].x;
					  image_disp[iloop].y= molecule[start_mol].y- molecule[start_mol+iloop].y;
					  image_disp[iloop].z= molecule[start_mol].z- molecule[start_mol+iloop].z;
				
					  min_image(&(image_disp[iloop].x),
						    &(image_disp[iloop].y), 
						    &(image_disp[iloop].z)); 

					  molecule[start_mol+iloop].x = molecule[start_mol].x - image_disp[iloop].x;
					  molecule[start_mol+iloop].y = molecule[start_mol].y - image_disp[iloop].y;
					  molecule[start_mol+iloop].z = molecule[start_mol].z - image_disp[iloop].z;
					}
				   }
*/

	/*********************************************/
	/*** Write frame to pdb file             *****/
	/*** copy atoms into pdb_frame and       *****/
	/*** only send on final pass to allow    *****/
	/*** correct placement of connect info   *****/
	/*********************************************/
                             if (anal_flags[ianal].pdb || anal_flags[ianal].msd)
                               {
                                  if (frame_index == 0 && this_molecule == 0 && !done_pdb)
                                    {
	                               if ( anal_flags[ianal].pdb && (pdb_fp = fopen (&pdb_file[0],"w") ) == NULL)
                                        {
                                         printf ("ERROR Opening pdb file %s for writing\n", pdb_file);
                                         perror ("open");
                                         return 1;
	                                }
                                      printf("Openned pdb file here\n");
                                      done_pdb=TRUE;
                                      num_pdb_atoms=0;
                                      num_pdb_writes=0;

                                      for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
				       {
        /*********************************************/
        /* In first frame move all atoms to min_image*/
        /* with first in list                        */
        /*********************************************/

        /* Get vector from current atom to first in list */

                                         vec[0]= molecule[iloop].x - molecule[start_mol].x;
                                         vec[1]= molecule[iloop].y - molecule[start_mol].y;
                                         vec[2]= molecule[iloop].z - molecule[start_mol].z;

        /* Find minimum image of this atom from first in list */

                                         min_image(&vec[0], &vec[1], &vec[2]); 

        /* Remember initial position of this atom          */
             
                                         vec1[0]= molecule[iloop].x;
                                         vec1[1]= molecule[iloop].y;
                                         vec1[2]= molecule[iloop].z;

        /* Move the atom to its minimum image co-ordinates */

                                         molecule[iloop].x = molecule[start_mol].x+vec[0];
                                         molecule[iloop].y = molecule[start_mol].y+vec[1];
                                         molecule[iloop].z = molecule[start_mol].z+vec[2];

        /* Record the displacement used for this atom      */
                                           
                                         image_disp[iloop].x = molecule[iloop].x - vec1[0];
                                         image_disp[iloop].y = molecule[iloop].y - vec1[1];
                                         image_disp[iloop].z = molecule[iloop].z - vec1[2];

                                         last_frame[iloop] = molecule[iloop];
                                       }
                                    }

	/*********************************************/
	/*** Remove periodic boundary for all ********/
	/*********************************************/
                                  if ( frame_index > 0 && !want_raw )
                                    {

        /* Work out displacement of this atom from last time step */
        /* Effective current position of atom is current co-ords  */
        /* plus image_disp vector                                 */

                                         dx = last_frame[start_mol].x - molecule[start_mol].x - image_disp[start_mol].x;
                                         dy = last_frame[start_mol].y - molecule[start_mol].y - image_disp[start_mol].y;
                                         dz = last_frame[start_mol].z - molecule[start_mol].z - image_disp[start_mol].z;

        /* Check for lattice vector size jumps that indicate change of min-image */

                                         if( dx > 0.5*latt_vec[0])
                                           {
                                             image_disp[start_mol].x += latt_vec[0];
                                             image_disp[start_mol].y += latt_vec[1];
                                             image_disp[start_mol].z += latt_vec[2];
                                           }

                                         if( dx <-0.5*latt_vec[0])
                                           {
                                             image_disp[start_mol].x -= latt_vec[0];
                                             image_disp[start_mol].y -= latt_vec[1];
                                             image_disp[start_mol].z -= latt_vec[2];
                                           }

                                         if( dy > 0.5*latt_vec[4])
                                           {
                                             image_disp[start_mol].x += latt_vec[3];
                                             image_disp[start_mol].y += latt_vec[4];
                                             image_disp[start_mol].z += latt_vec[5];
                                           }

                                         if( dy < -0.5*latt_vec[4])
                                           {
                                             image_disp[start_mol].x -= latt_vec[3];
                                             image_disp[start_mol].y -= latt_vec[4];
                                             image_disp[start_mol].z -= latt_vec[5];
                                           }

                                         if( dz > 0.5*latt_vec[8])
                                           {
                                             image_disp[start_mol].x += latt_vec[6];
                                             image_disp[start_mol].y += latt_vec[7];
                                             image_disp[start_mol].z += latt_vec[8];
                                           }

                                         if( dz < -0.5*latt_vec[8])
                                           {
                                             image_disp[start_mol].x -= latt_vec[6];
                                             image_disp[start_mol].y -= latt_vec[7];
                                             image_disp[start_mol].z -= latt_vec[8];
                                           }

                                          molecule[start_mol].x += image_disp[start_mol].x;
                                          molecule[start_mol].y += image_disp[start_mol].y;
                                          molecule[start_mol].z += image_disp[start_mol].z;

                                          last_frame[start_mol] = molecule[start_mol];

                                      for (iloop=start_mol+1; iloop <= demarcation[index].end; iloop++)
			               {
        /*********************************************/
        /* Move all other atoms to min_image         */
        /* with first in list                        */
        /*********************************************/

        /* Get vector from current atom to first in list */

                                         vec[0]= molecule[iloop].x - molecule[start_mol].x;
                                         vec[1]= molecule[iloop].y - molecule[start_mol].y;
                                         vec[2]= molecule[iloop].z - molecule[start_mol].z;

        /* Find minimum image of this atom from first in list */

                                         min_image(&vec[0], &vec[1], &vec[2]); 

        /* Remember initial position of this atom          */
             
                                         vec1[0]= molecule[iloop].x;
                                         vec1[1]= molecule[iloop].y;
                                         vec1[2]= molecule[iloop].z;

        /* Move the atom to its minimum image co-ordinates */

                                         molecule[iloop].x = molecule[start_mol].x+vec[0];
                                         molecule[iloop].y = molecule[start_mol].y+vec[1];
                                         molecule[iloop].z = molecule[start_mol].z+vec[2];

        /* Record the displacement used for this atom      */
                                           
                                         image_disp[iloop].x = molecule[iloop].x - vec1[0];
                                         image_disp[iloop].y = molecule[iloop].y - vec1[1];
                                         image_disp[iloop].z = molecule[iloop].z - vec1[2];

                                         last_frame[iloop] = molecule[iloop];

                                         }
        /*********************************************/
        /*** Work out centre of mass at new position */
        /*** for MSD analysis                        */
        /*********************************************/
			         if (anal_flags[ianal].msd && ( start_mol > 0 || num_molecules == 1 ) )
			           {

                                      if ( want_raw)
                                        {
                                           printf("MSD analysis of raw data is meaningless\n");
                                           exit(0);
                                        }
			              centre_of_mass(&c_of_m[0], &total_mass, &molecule[start_mol], 
				    			                    demarcation[index].num-1 );

                                      num_cofm_list++;
                                      cofm_list[num_cofm_list].x = c_of_m[0];
                                      cofm_list[num_cofm_list].y = c_of_m[1];
                                      cofm_list[num_cofm_list].z = c_of_m[2];

/*
                                      printf("cofm %10.6f %10.6f %10.6f\n", cofm_list[num_cofm_list].x,
                                                                            cofm_list[num_cofm_list].y,
                                                                            cofm_list[num_cofm_list].z);
*/
                                   }

                                 }
	/*********************************************/

                                  for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
				    {
                                      pdb_frame[num_pdb_atoms] = molecule[iloop];
                                      num_pdb_atoms++;
                                    }
                                  num_pdb_writes++;

                                  if (num_pdb_writes == tot_types)
                                    {
                                       write_pdb( pdb_fp, &pdb_frame[0], &abc[0], num_pdb_atoms);
                                       num_pdb_atoms=0;
                                       num_pdb_writes=0;
                                    }
                               }
	/*********************************************/
	/*** Insert additional atoms             *****/
	/*********************************************/

			    if (anal_flags[ianal].add_atom)
			      {
	/********************************************/
	/*** Crude implementation for TIP4P case ****/
	/********************************************/
				bis_found = 0;
				bis_index1 = -1;
				bis_index2 = -1;
				bis_index3 = -1;
				for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
				  {
				    if ( strcmp(molecule[iloop].label, bisector_atoms[0].label ) == 0)
				       { 
					  if   (bis_index1 < 0)  
					     {
						bis_index1 = iloop;
						bis_found++;
					     }
					  else if (bis_index3 < 0) 
					     {
						bis_index3 = iloop;
						bis_found++;
					     }
				       }
				    if ( strcmp(molecule[iloop].label, bisector_atoms[1].label ) == 0)
				       {
					  bis_index2 = iloop;
					  bis_found++; 
				       }
				    if (bis_found == 3) break;
				  }

				if (bis_found == 3) 
				  {
				     printf("Found bisector triplet %d (%s) %d (%s) %d (%s)\n", 
						   bis_index1, molecule[bis_index1].label, 
						   bis_index2, molecule[bis_index2].label,
						   bis_index3, molecule[bis_index3].label);

	/*********************************************/
	/*** Adjust for TIP4P*************************/
	/*********************************************/
		       
				       oh_dist = 0.9572;
				       oh_angl = 104.52/RAD_TO_DEG;
	    
	/* bond vectors */
	    
				       vec[0] = molecule[bis_index1].x - molecule[bis_index2].x;
				       vec[1] = molecule[bis_index1].y - molecule[bis_index2].y;
				       vec[2] = molecule[bis_index1].z - molecule[bis_index2].z;

				       unit_vector(&vec[0]);

				       vec1[0] = molecule[bis_index3].x - molecule[bis_index2].x;
				       vec1[1] = molecule[bis_index3].y - molecule[bis_index2].y;
				       vec1[2] = molecule[bis_index3].z - molecule[bis_index2].z;

				       unit_vector(&vec1[0]);

				       vec_cross(&vec[0], &vec1[0], &norm[0]);
	    
				       unit_vector(&norm[0]);
	    
				       vec_cross(&norm[0], &vec[0], &perp[0]);
	    
				       unit_vector(&perp[0]);
	    
				       vec1[0] = vec[0]*cos(oh_angl) + perp[0]*sin(oh_angl);
				       vec1[1] = vec[1]*cos(oh_angl) + perp[1]*sin(oh_angl);
				       vec1[2] = vec[2]*cos(oh_angl) + perp[2]*sin(oh_angl);
	    
	/* Set angle right */
	    
				       molecule[bis_index1].x = molecule[bis_index2].x + oh_dist * vec[0];
				       molecule[bis_index1].y = molecule[bis_index2].y + oh_dist * vec[1];
				       molecule[bis_index1].z = molecule[bis_index2].z + oh_dist * vec[2];
	    
				       molecule[bis_index3].x = molecule[bis_index2].x + oh_dist * vec1[0];
				       molecule[bis_index3].y = molecule[bis_index2].y + oh_dist * vec1[1];
				       molecule[bis_index3].z = molecule[bis_index2].z + oh_dist * vec1[2];
	    
	/*********************************************************/
	/** Work out co-ordinates of new atom along the bisector */
	/*********************************************************/

	/* bisector is along sum of unit vector directions */

				     vec[0] += vec1[0];
				     vec[1] += vec1[1];
				     vec[2] += vec1[2];

				     unit_vector(&vec[0]);

				     vec[0] = bisect_disp*vec[0];
				     vec[1] = bisect_disp*vec[1];
				     vec[2] = bisect_disp*vec[2];


	/* Define new atom */
				     strcpy(new_atom.label, "M");
				     new_atom.x = molecule[bis_index2].x + vec[0];
				     new_atom.y = molecule[bis_index2].y + vec[1];
				     new_atom.z = molecule[bis_index2].z + vec[2];

	/* Update atom list */
				     (demarcation[index].end)++;
				     (demarcation[index].num)++;

				     max_atom_index++;
			for (iloop=max_atom_index; iloop > demarcation[index].end; iloop--) molecule[iloop] = molecule[iloop-1];                                  
				     printf("Inserting M at list position %d\n", demarcation[index].end);
				     molecule[demarcation[index].end] = new_atom; 

				     for (iloop = index+1; iloop <= tot_types; iloop++)
				       {
					  (demarcation[iloop].start)++;
					  (demarcation[iloop].end)++;
					  (demarcation[iloop].num)++;
				       }
				  }
				else
				  {
				    printf("Could not find a triplet for add_atom in this case\n");
				  }
			      }

			    if (anal_flags[ianal].coords)
			      {
			    printf("Frame %d, molecule %d starts at %d ends at %d\n", 
                                        frame_index, index, start_mol, demarcation[index].end);

	/*********************************************/
	/* HAVE FLAG TO GIVE CO_ORDS *****************/
	/* DO NOT GO TO COFM IF NO MASSES SET  *******/
	/*********************************************/

        /*** Just output all molecules when we get to the first one for complete frames in arc file ***/
                            if (this_molecule == 0 && itype == 0)
                               {
/*
                           if (itype==1) 
                            {
                              printf("max_atom_index=%d\n",max_atom_index);
                              exit(0);
                            }
*/
                                 if (frame_index == 0) 
                                   {
                                      start_frame=TRUE;
                                   
                                 sprintf(car_file,"movie.arc");
                                 sprintf(c_title_line,"frame number %d from dlpoly run\n",frame_index+1);

        /** flag header, date, and title as blank so defaults are used ***/
                                 header_line[0]=-1;
                                 title_line[0]=-1;
                                 date_line[0]=-1;

        /** Assume 1 by 1 by 1 cell should be output ****/
                                 super[0] = 1;
                                 super[1] = 1;
                                 super[2] = 1;

	                           if ( (car_fp = fopen (&car_file[0],"w") ) == NULL)
                                     {
                                        printf("ERROR openning car file for writing\n");
                                        exit(0);
                                     }
                                   }
                                 else 
                                     {
                                        start_frame=FALSE;
                                     }

                                 for (iloop=0; iloop < max_atom_index; iloop++)
                                   {
                                       mol_nums[iloop] = 0;
                                       fix_flags[iloop].fx=FALSE;
                                       fix_flags[iloop].fy=FALSE;
                                       fix_flags[iloop].fz=FALSE;
                                       molecule[iloop].mol=0;
                                   }

                                 gather_molecule(&molecule[start_mol], max_atom_index, 0);

                                 write_car( car_fp, &header_line[0], &title_line[0], &c_title_line[0],
		                            &date_line[0], &molecule[start_mol], &mol_nums[0], 
                                            pbc, &abc[0], max_atom_index+1, 1.0, 
                                            start_frame, &super[0], &latt_vec[0], &recip_latt_vec[0], 
                                            &fix_flags[0]);

				 for (iloop=start_mol; iloop <= max_atom_index; iloop++)
				    {
				       printf("%3d %8s %10.6f %10.6f %10.6f\n", iloop, molecule[iloop].label,
					     molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);
				    }
                                 }
			       }

	/*********************************************/
	/** Write out in hin file format *************/
	/*********************************************/

			    if (anal_flags[ianal].hin)
			      {
                              }

	/*********************************************/
	/* Calculate close contacts with all         */
	/* other molecules if requested              */
	/*********************************************/

			   if (anal_flags[ianal].contact)
			     {
			      j_index=0;
			      for (jtype=0; jtype < num_molecules; jtype++)
				 {
				  for (j_this_mol=0; j_this_mol < num_this_mol[jtype]; j_this_mol++)
				    {
	/********************************************/
	/* Do not look for intramolecular contacts  */
	/********************************************/

				       if ( index != j_index )
					 {
					   j_start_mol= demarcation[j_index].start;   
	  
					   close_contact(&molecule[start_mol], 
							 demarcation[index].num-1,
							 &molecule[j_start_mol], 
							 demarcation[j_index].num-1,
							 cutoff2, pbc, &num_cc,
							 &clash_list[0]);

					   if (num_cc >= 0)
					     {
						printf("\nFound %d close contact for molecules %d and %d\n",
							  num_cc+1, index+1, j_index+1);
						for (iloop=0; iloop <= num_cc; iloop++)
						  {
						    printf("Atom %4d (%4s) is %10.4fA from atom %4d (%4s)\n",
								  1+start_mol+clash_list[iloop].index1,
								  clash_list[iloop].label1,
								  clash_list[iloop].sep,
								  1+j_start_mol+clash_list[iloop].index2,
								  clash_list[iloop].label2);
							    
						  }
					     }
					 }
				       j_index++;
				     }
				 }
			      }

	/*********************************************/
	/* Work out dipole moment for each molecule  */
	/*********************************************/

			    if (anal_flags[ianal].dipoles )
			      {
				mol_dipole( &molecule[start_mol], 
					    demarcation[index].num-1, 
					    &mole_dips[index] );

				printf("Molecular dipole : %10.6f %10.6f %10.6f",
					  mole_dips[index].x, mole_dips[index].y, mole_dips[index].z); 

				dip_size = mole_dips[index].x*mole_dips[index].x 
					 + mole_dips[index].y*mole_dips[index].y
					 + mole_dips[index].z*mole_dips[index].z;

				dip_size = sqrt(dip_size);

				printf("  magnitude: %10.6f\n", dip_size); 
			      }

	/********************************************/
	/*** Work out forster factor for each pair **/
	/********************************************/
			   if (anal_flags[ianal].forster )
			     {
			      j_index=0;

			      for (jtype=0; jtype < num_molecules; jtype++)
				 {
				  for (j_this_mol=0; j_this_mol < num_this_mol[jtype]; j_this_mol++)
				    {
	/***********************************************/
	/* Do not look for intramolecular interactions */
	/***********************************************/

				       if ( use_type[jtype] && index < j_index )
					 {
					   j_start_mol= demarcation[j_index].start;
	 
					   forster_kappa(&molecule[start_mol],
							 demarcation[index].num-1,
							 &molecule[j_start_mol],
							 demarcation[j_index].num-1,
							 &kappas[num_kappas]);

					   tot_kappa += kappas[num_kappas]*kappas[num_kappas];

			  printf("Forster kappa for molecule %d with %d = %10.6f kappa squared = %10.6f\n",
								    index, j_index, kappas[num_kappas], 
								    kappas[num_kappas]* kappas[num_kappas]);
					   num_kappas++;
					 }
				       j_index++;
				     }
				 }
			      }

			    if (anal_flags[ianal].mo_inertia)
			      {

	/*********************************************/
	/* Work out Centre of Mass for each molecule */
	/*********************************************/

			    centre_of_mass(&c_of_m[0], &total_mass, &molecule[start_mol], 
							       demarcation[index].num-1 );
	/******************************/
	/* Work out Moment of Inertia */
	/******************************/

			    moments_of_inertia(&molecule[start_mol], 
					       demarcation[index].num-1, &c_of_m[0],
					       &m_of_inertia[0], &moi_eigenvals[0] );

			    big_eigen=0.0;
			    for (iloop=0; iloop < 3; iloop++)
			      {
				eigen_vec_3b3( &m_of_inertia[0], moi_eigenvals[iloop], 
							     &moi_eigenvecs[3*iloop] );

			      }

	/***************************************/
	/* Work out order of moi eigenvalues ***/
	/***************************************/

			     big_eigen= 0.0;
			     small_eigen= 1E28;
			     for (iloop=0; iloop < 3; iloop++)
			      {
				if (moi_eigenvals[iloop] > big_eigen) 
				  {
				    big_eigen=moi_eigenvals[iloop];
				    big_index= iloop;
				  }
				if (moi_eigenvals[iloop] < small_eigen)
				  {
				    small_eigen= moi_eigenvals[iloop];
				    small_index= iloop;
				  }
			      }

			     map_moi_order[0]=big_index;
			     map_moi_order[2]=small_index;
			  
			     for (iloop=0; iloop < 3; iloop++)
				   if (iloop != big_index && iloop != small_index) 
							    map_moi_order[1] = iloop;

	/*****************************************************/
	/*** Print out components required *******************/
	/*****************************************************/
	/*****************************************************/
	/*** Deal with smallest moi eigenvalue ***************/
	/*****************************************************/
			  
			    if (anal_flags[ianal].big_med_small[0])
			      {
				index_new= 3*map_moi_order[0];

				if (anal_flags[ianal].angles)
				  {
				     angle[0]= acos(moi_eigenvecs[index_new]);
				     angle[1]= acos(moi_eigenvecs[index_new+1]);
				     angle[2]= acos(moi_eigenvecs[index_new+2]);

				     if (angle[2] > 90.0) angle[2] = 180.0 - angle[2];

				     printf("Largest eigenvalue of moi = %10.6f : eigenvec: %10.6f %10.6f %10.6f angles with xyz: %10.6f %10.6f %10.6f\n",
					      moi_eigenvals[map_moi_order[0]], moi_eigenvecs[index_new],
					      moi_eigenvecs[index_new+1], moi_eigenvecs[index_new+2],
					      angle[0], angle[1], angle[2]);
				  }
				else
				  {
				     printf("Largest eigenvalue of moi = %10.6f : eigenvec: %10.6f %10.6f %10.6f\n",
					      moi_eigenvals[map_moi_order[0]], moi_eigenvecs[index_new],
					      moi_eigenvecs[index_new+1], moi_eigenvecs[index_new+2]);
				  }

			       if (anal_flags[ianal].sum_comps)
				  {
				     all_mols_moi[0].totx += moi_eigenvecs[index_new]*moi_eigenvecs[index_new];
				     all_mols_moi[0].toty += moi_eigenvecs[index_new+1]*moi_eigenvecs[index_new+1];
				     all_mols_moi[0].totz += moi_eigenvecs[index_new+2]*moi_eigenvecs[index_new+2];
				     (all_mols_moi[0].num)++;
				  }
	/********************************************/
	/** Store moi eigenvector for correlation ***/
	/********************************************/
			       if (anal_flags[ianal].correlation)
				  {
				     la_moi_vs_time[la_icorr].x = moi_eigenvecs[index_new];
				     la_moi_vs_time[la_icorr].y = moi_eigenvecs[index_new+1];
				     la_moi_vs_time[la_icorr].z = moi_eigenvecs[index_new+2];
				     la_icorr++;
				  }
			      }

	/*****************************************************/
	/*** Deal with medium moi eigenvalue *****************/
	/*****************************************************/
			    if (anal_flags[ianal].big_med_small[1])
			      {
				index_new= 3*map_moi_order[1];
	 
				if (anal_flags[ianal].angles)
				  {
				    angle[0]= acos(moi_eigenvecs[index_new]);
				    angle[1]= acos(moi_eigenvecs[index_new+1]);
				    angle[2]= acos(moi_eigenvecs[index_new+2]);

				     if (angle[2] > 90.0) angle[2] = 180.0 - angle[2];

				    printf("Middle eigenvalue of moi = %10.6f : eigenvec: %10.6f %10.6f %10.6f angles with xyz: %10.6f %10.6f %10.6f\n",
					     moi_eigenvals[map_moi_order[1]], moi_eigenvecs[index_new],
					     moi_eigenvecs[index_new+1], moi_eigenvecs[index_new+2],
					     angle[0], angle[1], angle[2]);
				  }
				else
				  {
				    printf("Middle eigenvalue of moi = %10.6f : eigenvec: %10.6f %10.6f %10.6f\n",
					     moi_eigenvals[map_moi_order[1]], moi_eigenvecs[index_new],
					     moi_eigenvecs[index_new+1], moi_eigenvecs[index_new+2]);
				  }
			       if (anal_flags[ianal].sum_comps)
				  {
				     all_mols_moi[1].totx += moi_eigenvecs[index_new]*moi_eigenvecs[index_new];
				     all_mols_moi[1].toty += moi_eigenvecs[index_new+1]*moi_eigenvecs[index_new+1];
				     all_mols_moi[1].totz += moi_eigenvecs[index_new+2]*moi_eigenvecs[index_new+2];
				     (all_mols_moi[1].num)++;
				  }
	/********************************************/
	/** Store moi eigenvector for correlation ***/
	/********************************************/
			       if (anal_flags[ianal].correlation)
				  {
				     me_moi_vs_time[me_icorr].x = moi_eigenvecs[index_new];
				     me_moi_vs_time[me_icorr].y = moi_eigenvecs[index_new+1];
				     me_moi_vs_time[me_icorr].z = moi_eigenvecs[index_new+2];
				     me_icorr++;
				  }

			      }

	/*****************************************************/
	/*** Deal with biggest moi eigenvalue ****************/
	/*****************************************************/
			    if (anal_flags[ianal].big_med_small[2])
			      {
				index_new= 3*map_moi_order[2];

				if (anal_flags[ianal].angles)
				  { 
				    angle[0]= acos(moi_eigenvecs[index_new]);
				    angle[1]= acos(moi_eigenvecs[index_new+1]);
				    angle[2]= acos(moi_eigenvecs[index_new+2]);

				     if (angle[2] > 90.0) angle[2] = 180.0 - angle[2];

				    printf("Small eigenvalue of moi = %10.6f : eigenvec: %10.6f %10.6f %10.6f angles with xyz: %10.6f %10.6f %10.6f\n",
					    moi_eigenvals[map_moi_order[2]], moi_eigenvecs[index_new],
					    moi_eigenvecs[index_new+1], moi_eigenvecs[index_new+2],
					    angle[0], angle[1], angle[2]);
				  }
				else
				  {
				    printf("Small eigenvalue of moi = %10.6f : eigenvec: %10.6f %10.6f %10.6f\n",
					    moi_eigenvals[map_moi_order[2]], moi_eigenvecs[index_new],
					    moi_eigenvecs[index_new+1], moi_eigenvecs[index_new+2]);
				  }
			       if (anal_flags[ianal].sum_comps)
				  {
				     all_mols_moi[2].totx += moi_eigenvecs[index_new]*moi_eigenvecs[index_new];
				     all_mols_moi[2].toty += moi_eigenvecs[index_new+1]*moi_eigenvecs[index_new+1];
				     all_mols_moi[2].totz += moi_eigenvecs[index_new+2]*moi_eigenvecs[index_new+2];
				     (all_mols_moi[2].num)++;
				  }
	/********************************************/
	/** Store moi eigenvector for correlation ***/
	/********************************************/
			       if (anal_flags[ianal].correlation)
				  {
				     sm_moi_vs_time[sm_icorr].x = moi_eigenvecs[index_new];
				     sm_moi_vs_time[sm_icorr].y = moi_eigenvecs[index_new+1];
				     sm_moi_vs_time[sm_icorr].z = moi_eigenvecs[index_new+2];
				     sm_icorr++;
				  }
			      }
			 }

          /*******************************************************************/
          /*** Velocity auto-correlation data ********************************/
          /*** Note that currently the correlation flag is shared with moi ***/
          /*******************************************************************/

                      if (anal_flags[ianal].correlation)
                         {

                            printf("Doing velocity list for molecule starting at %d and ending at %d\n", 
                                                                  start_mol, demarcation[index].end );
                            vel_list[vel_icorr].x = 0.0;
                            vel_list[vel_icorr].y = 0.0;
                            vel_list[vel_icorr].z = 0.0;

                            for ( iatom = start_mol; iatom <= demarcation[index].end; iatom++ )
                               {
                                 vel_list[vel_icorr].x += molecule[iatom].vx;
                                 vel_list[vel_icorr].y += molecule[iatom].vy;
                                 vel_list[vel_icorr].z += molecule[iatom].vz;
                               }
/****************************************/
/*** Molecular velocity is average of ***/
/*** members velocities               ***/
/****************************************/
                             vel_list[vel_icorr].x = vel_list[vel_icorr].x / (double) demarcation[index].num;
                             vel_list[vel_icorr].y = vel_list[vel_icorr].y / (double) demarcation[index].num;
                             vel_list[vel_icorr].z = vel_list[vel_icorr].z / (double) demarcation[index].num;
 
/*
                            printf("\n\nTotal velocity of molecule = %10.6f, %10.6f, %10.6f.\n", vel_list[vel_icorr].x,
                                                                                                 vel_list[vel_icorr].y,
                                                                                                 vel_list[vel_icorr].z);
*/

                            vel_icorr++;
                         }

/*** Analysis added for co-crystals *****/

                      if (anal_flags[ianal].bonds)
                         {


                         for ( iatom = start_mol; iatom <= demarcation[index].end; iatom++ )
                           {

                        do_this = FALSE;
                        for ( iii = 0; iii <= num_std_bnds; iii++) 
                           if (strcmp(molecule[iatom].elem, monit_bnds[iii].name)==0)
                                 do_this = TRUE;
                     if (do_this) 
                       {
                        for (ineigh=0; ineigh <= molecule[iatom].num_neigh; ineigh++) 
                           {
                               this_neigh = molecule[iatom].neighb[ineigh];

                        do_this = FALSE;
                        for ( iii = 0; iii <= num_std_bnds; iii++) 
                           if (strcmp(molecule[this_neigh].elem, monit_bnds[iii].name)==0)
                                 do_this = TRUE;

                               if (iatom > this_neigh && do_this)
                                 {
                          /*    printf("measuring %d %s to %d %s = ",               */
                          /*                iatom, molecule[iatom].label,           */
                          /*                this_neigh, molecule[this_neigh].label);*/

                                   dx =  molecule[this_neigh].x - molecule[iatom].x;
                                   dy =  molecule[this_neigh].y - molecule[iatom].y;
                                   dz =  molecule[this_neigh].z - molecule[iatom].z;

                                   min_image(&dx, &dy, &dz);

                                   dist = sqrt(dx*dx +dy*dy +dz*dz); 

                           /*      printf("%10.6f\n",dist); */

/*** Work out which bin this belongs to ****/

                                   guess = (dist-bnd_low)/bnd_delta - 0.5;
 
                                   iguess = guess;
                                   iguess++;

                           /*      printf("iguess = %d, low %10.6f high %10.6f\n",         */
                           /*                       iguess, (iguess-0.5)*bnd_delta+bnd_low,*/
                           /*                       (iguess+0.5)*bnd_delta+bnd_low );      */

/**** Increment that bins contents *********/

                                   bins[iguess]++;
                                 }
                             }
                            }
                           }

                         } 

		      index++;
		     }
		   }
		else
		   {
	/*****************************************************/
	/*** Even if not using need to keep track of index ***/
	/*****************************************************/
		      for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
											     index++;
		   }
		}
	/*********************************************/
	/*** Report frame averages *******************/
	/*********************************************/
	       if (anal_flags[ianal].forster )
		{
		  printf("Average kappa squared for frame = %10.6f\n",
					   tot_kappa/num_kappas);
		}

	/*********************************************/
	/*** Report sum and average moi components ***/
	/*********************************************/
	      if (anal_flags[ianal].mo_inertia)
		{
		  if (anal_flags[ianal].sum_comps)
		    {
		       if (anal_flags[ianal].big_med_small[0])
			 {
			    printf("Total square for all components of eigenvector corresponding to Largest moi eigenvector\n");
			    printf("Total x : %10.6f, Total y : %10.6f, Total z : %10.6f\n",
				      all_mols_moi[0].totx, all_mols_moi[0].toty, all_mols_moi[0].totz);
			    printf("Totals taken over %d molecule/frame combinations\n", all_mols_moi[0].num);
			    rnum= all_mols_moi[0].num;
			    printf("Root mean square x: %10.6f, y: %10.6f, z: %10.6f\n\n",
				      sqrt(all_mols_moi[0].totx/rnum), sqrt(all_mols_moi[0].toty/rnum), sqrt(all_mols_moi[0].totz/rnum));
			 }
		       if (anal_flags[ianal].big_med_small[1])
			 {
			    printf("Totals square for all components of eigenvector corresponding to Medium moi eigenvector\n");
			    printf("Total x : %10.6f, Total y : %10.6f, Total z : %10.6f\n",
				      all_mols_moi[1].totx, all_mols_moi[1].toty, all_mols_moi[1].totz);
			    printf("Totals taken over %d molecule/frame combinations\n", all_mols_moi[1].num);
			    rnum= all_mols_moi[1].num;
			    printf("Root mean square x: %10.6f, y: %10.6f, z: %10.6f\n\n",
				      sqrt(all_mols_moi[1].totx/rnum), sqrt(all_mols_moi[1].toty/rnum), sqrt(all_mols_moi[1].totz/rnum));
			 }
		       if (anal_flags[ianal].big_med_small[2])
			{
			    printf("Totals square for all components of eigenvector corresponding to Smallest moi eigenvector\n");
			    printf("Total x : %10.6f, Total y : %10.6f, Total z : %10.6f\n",
				      all_mols_moi[2].totx, all_mols_moi[2].toty, all_mols_moi[2].totz);
			    printf("Totals taken over %d molecule/frame combinations\n", all_mols_moi[2].num);
			    rnum= all_mols_moi[2].num;
			    printf("Root mean square x: %10.6f, y: %10.6f, z: %10.6f\n\n",
				      sqrt(all_mols_moi[2].totx/rnum), sqrt(all_mols_moi[2].toty/rnum), sqrt(all_mols_moi[2].totz/rnum));
			}
		   }
		}
	/*********************************************/
	/*** Print out adapted atom list *************/
	/*********************************************/

	if (anal_flags[ianal].add_atom)
	    {
	       printf("New atom list:\n");
	       for ( iatom=0; iatom <= max_atom_index; iatom++)
		 {
		    printf("%-8s %d\n", molecule[iatom].label, iatom+1);
		    printf("%20.8f%20.8f%20.8f\n",
				    molecule[iatom].x, molecule[iatom].y, molecule[iatom].z);
		 }
	    }

	    }

          if (anal_flags[ianal].coords) fclose(car_fp);

	  if (have_config) break;
	 }

	/***********************/
	/* Closes HISTORY File */
	/***********************/

	if ( (fclose (history_fp) ) == EOF)
	    {
	    printf ("%s ERROR Closing Input File\n");
	    perror ("close");
	    return 3;
	    }

	/**********************************/
	/*** report any histograms ********/
	/**********************************/

       centre = bnd_low;
       for (ibin = 0; ibin < num_bins; ibin++) 
        {
          printf("%d, %10.6f, %d\n",ibin,centre,bins[ibin]);
          centre += bnd_delta;
       }

	/**********************************/
	/** Do any correlations required **/
	/**********************************/

	for (ianal=0; ianal <= num_anal_flags; ianal++)
	   {
/************************************************/
/*** Calculate MSD plot *************************/
/************************************************/
             if (anal_flags[ianal].msd)
               {
                 if ( num_molecules > 1)
                   {
                     msd_calc(&cofm_list[0], &msd[0], &time_dots[0],
                              time_step, min_msd_ind_diff, 
                              max_msd_ind_diff, frame_index, num_this_mol[1]-1 );

                     printf("Back from msd_calc\n");
                   }
                 else
                   {
                     msd_calc(&cofm_list[0], &msd[0], &time_dots[0],
                              time_step, min_msd_ind_diff, 
                              max_msd_ind_diff, frame_index, num_this_mol[0]-1 );
                   }

                  if ( (msd_fp = fopen ("msd.csv","w") ) == NULL)
                    {
                      printf ("ERROR Opening msd.csv FILE\n");
                      return 1;
                    }
/*
                  for ( iloop=0; iloop<= max_msd_ind_diff-min_msd_ind_diff; iloop++)
                    {
                      printf("msd time %10.6f  value %10.6f\n", time_dots[iloop], msd[iloop]);
                    }
*/
                  write_csv(msd_fp, "time", "msd", "msd", 
                                 &time_dots[0] , &msd[0], &msd[0], FALSE, 
                                max_msd_ind_diff-min_msd_ind_diff);

                  printf("Written msd csv file closing file pointer\n");
                 
                  fclose(msd_fp);
               }

	     if (anal_flags[ianal].big_med_small[0]) 
	       {
                 if ( num_molecules > 1)
                   {
		     dot_correlation(&la_moi_vs_time[0], &la_moi_dots[0],
                             &time_dots[0], time_step, min_corr_ind_diff, 
                             max_corr_ind_diff,
                             frame_index, num_this_mol[1]-1 ); 
                   }
                 else
                   {
		     dot_correlation(&la_moi_vs_time[0], &la_moi_dots[0],
                             &time_dots[0], time_step, min_corr_ind_diff, 
                             max_corr_ind_diff,
                             frame_index, num_this_mol[0]-1 ); 
                   }

                  printf("Dot product correlation for largest eigenvalue of moi\n");
                  index=0;
                  for (iloop= min_corr_ind_diff; iloop <= max_corr_ind_diff; iloop++)
                     {
                       printf(" %d   %10.6f\n", iloop, la_moi_dots[index]);
                       index++;
                     } 
               }
   
             if (anal_flags[ianal].big_med_small[1])
               {
                  dot_correlation(&me_moi_vs_time[0], &me_moi_dots[0],
                                  &time_dots[0], time_step,
                                  min_corr_ind_diff, max_corr_ind_diff,
                                  frame_index, num_this_mol[1]-1 );

                  printf("Dot product correlation for medium eigenvalue of moi\n");
                  index=0;
                  for (iloop= min_corr_ind_diff; iloop <= max_corr_ind_diff; iloop++)
                     {
                       printf(" %d   %10.6f\n", iloop, me_moi_dots[index]);
                       index++;
                     } 
               }

             if (anal_flags[ianal].big_med_small[2])
               {
                  dot_correlation(&sm_moi_vs_time[0], &sm_moi_dots[0],
                                  &time_dots[0], time_step,
                                  min_corr_ind_diff, max_corr_ind_diff,
                                  frame_index, num_this_mol[1]-1 );

                  printf("Dot product correlation for smallest eigenvalue of moi\n");
                  index=0;
                  for (iloop= min_corr_ind_diff; iloop <= max_corr_ind_diff; iloop++)
                     {
                       printf(" %d   %10.6f\n", iloop, sm_moi_dots[index]);
                       index++;
                     } 
               }

             if (anal_flags[ianal].correlation)
               {
                 if ( num_molecules > 1)
                   {
                  dot_correlation(&vel_list[0], &vel_dots[0],
                                  &time_dots[0], time_step,
                                  min_corr_ind_diff, max_corr_ind_diff,
                                  frame_index, num_this_mol[1]-1 );
                   }
                 else
                   {
                  dot_correlation(&vel_list[0], &vel_dots[0],
                                  &time_dots[0], time_step,
                                  min_corr_ind_diff, max_corr_ind_diff,
                                  frame_index, num_this_mol[0]-1 );
                   }

                  if ( (correlation_fp = fopen ("vel_correlation.csv","w") ) == NULL)
                    {
                      printf ("ERROR Opening velocity correlation FILE\n");
                      return 1;
                    }

                  write_csv(correlation_fp, "time", "velocity correlation", "velocity correlation", 
                                 &time_dots[0] , &vel_dots[0], &vel_dots[0], FALSE, 
                                max_corr_ind_diff-min_corr_ind_diff);

                  fclose(correlation_fp);
/*
                  printf(" Velocity dot product correlation function\n");
                  index=0;
                  for (iloop= min_corr_ind_diff; iloop <= max_corr_ind_diff; iloop++)
                     {
                       printf(" %d   %10.6f\n", iloop, vel_dots[index]);
                       index++;
                     }
*/
               }
           }

/********************************/
/* Closes Specified Output File */
/********************************/

/* if ( (fclose (fptr) ) == EOF)                   */
/*     {                                           */
/*     printf ("%s ERROR Closing Output File\n");  */
/*     perror ("close");                           */
/*     return 4;                                   */
/*     }                                           */
return 0;
}
