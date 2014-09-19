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

void write_car( int *p_header_line, int *p_title_line, int *p_date_line,
                atom *p_molecule, int pbc, double *p_abc, int num_atoms,
                int do_header);

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
                double *p_timestep, int *p_have_rdf, int *p_have_out);

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
                 int *p_num_rdfs,
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
int icount;
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

FILE *field_fp;
FILE *rdf_fp;
FILE *pdb_fp;
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

int num_kappas;

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
vector vel_list[MAX_CORR_VECS];

int frame_index, sm_icorr,me_icorr,la_icorr;
int vel_icorr=0;
double sm_moi_dots[MAXFRAMES];
double me_moi_dots[MAXFRAMES];
double la_moi_dots[MAXFRAMES];
double vel_dots[MAXFRAMES];

double eng_tot[MAXFRAMES];
double enthalpy[MAXFRAMES],   stat_time[MAXFRAMES]; 
double temp_tot[MAXFRAMES],   eng_cfg[MAXFRAMES]; 
double eng_vdw[MAXFRAMES],    eng_cou[MAXFRAMES]; 
double eng_bnd[MAXFRAMES],    eng_ang[MAXFRAMES]; 
double eng_dih[MAXFRAMES],    eng_tet[MAXFRAMES]; 
double time[MAXFRAMES],       eng_pv[MAXFRAMES];   
double temp_rot[MAXFRAMES],   vir_cfg[MAXFRAMES];  
double vir_vdw[MAXFRAMES],    vir_cou[MAXFRAMES];   
double vir_bnd[MAXFRAMES],    vir_ang[MAXFRAMES];   
double vir_con[MAXFRAMES],    vir_tet[MAXFRAMES];   
double volume[MAXFRAMES],     temp_shl[MAXFRAMES];
double eng_shl[MAXFRAMES],    vir_shl[MAXFRAMES];
double alpha[MAXFRAMES],      beta[MAXFRAMES];
double gamma[MAXFRAMES],      vir_pmf[MAXFRAMES];
double press[MAXFRAMES];

int num_in_list, have_statis;
int num_stat_list;

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
                       &have_rdf, &have_out );

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
                &num_rdfs,
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

		  if (anal_flags[iloop].mo_inertia && anal_flags[iloop].sum_comps)
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

		  printf("Carrying out Analysis pass number %d\n", ianal+1);

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

	/*********************************************/
	/*** Consider each of this molecule type *****/
	/*********************************************/

			   for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
			      {
				start_mol= demarcation[index].start;

	/*********************************************/
	/*** Write frame to pdb file             *****/
	/*** copy atoms into pdb_frame and       *****/
	/*** only send on final pass to allow    *****/
	/*** correct placement of connect info   *****/
	/*********************************************/
                             if (anal_flags[ianal].pdb || anal_flags[ianal].msd)
                               {
                                  if (frame_index == 0)
                                    {
	                               if ( anal_flags[ianal].pdb && (pdb_fp = fopen (&pdb_file[0],"w") ) == NULL)
                                        {
                                         printf ("ERROR Opening pdb file %s for writing\n", pdb_file);
                                         perror ("open");
                                         return 1;
	                                }
                                      num_pdb_atoms=0;
                                      num_pdb_writes=0;

                                      for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
				       {
        /*********************************************/
        /* In first frame move all atoms to min_image*/
        /* with first in list                        */
        /*********************************************/
                                         image_disp[iloop].x = 0.0;
                                         image_disp[iloop].y = 0.0;
                                         image_disp[iloop].z = 0.0;

                                         last_frame[iloop] = molecule[iloop];
                                       }
                                    }

	/*********************************************/
	/*** Remove periodic boundary for all ********/
	/*********************************************/
                                  if ( frame_index > 0 )
                                    {
                                      for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
			               {
                                         dx = last_frame[iloop].x - molecule[iloop].x - image_disp[iloop].x;
                                         dy = last_frame[iloop].y - molecule[iloop].y - image_disp[iloop].y;
                                         dz = last_frame[iloop].z - molecule[iloop].z - image_disp[iloop].z;

                                         if( dx > 0.9*latt_vec[0])
                                           {
                                             image_disp[iloop].x += latt_vec[0];
                                             image_disp[iloop].y += latt_vec[1];
                                             image_disp[iloop].z += latt_vec[2];
                                           }

                                         if( dx <-0.9*latt_vec[0])
                                           {
                                             image_disp[iloop].x -= latt_vec[0];
                                             image_disp[iloop].y -= latt_vec[1];
                                             image_disp[iloop].z -= latt_vec[2];
                                           }

                                         if( dy > 0.9*latt_vec[4])
                                           {
                                             image_disp[iloop].x += latt_vec[3];
                                             image_disp[iloop].y += latt_vec[4];
                                             image_disp[iloop].z += latt_vec[5];
                                           }

                                         if( dy < -0.9*latt_vec[4])
                                           {
                                             image_disp[iloop].x -= latt_vec[3];
                                             image_disp[iloop].y -= latt_vec[4];
                                             image_disp[iloop].z -= latt_vec[5];
                                           }

                                         if( dz > 0.9*latt_vec[8])
                                           {
                                             image_disp[iloop].x += latt_vec[6];
                                             image_disp[iloop].y += latt_vec[7];
                                             image_disp[iloop].z += latt_vec[8];
                                           }

                                         if( dz < -0.9*latt_vec[8])
                                           {
                                             image_disp[iloop].x -= latt_vec[6];
                                             image_disp[iloop].y -= latt_vec[7];
                                             image_disp[iloop].z -= latt_vec[8];
                                           }

                                          molecule[iloop].x += image_disp[iloop].x;
                                          molecule[iloop].y += image_disp[iloop].y;
                                          molecule[iloop].z += image_disp[iloop].z;

                                          last_frame[iloop] = molecule[iloop];
                                         }
        /*********************************************/
        /*** Work out centre of mass at new position */
        /*** for MSD analysis                        */
        /*********************************************/
			         if (anal_flags[ianal].msd && ( start_mol > 0 || num_molecules == 1 ) )
			           {
			              centre_of_mass(&c_of_m[0], &total_mass, &molecule[start_mol], 
				    			                    demarcation[index].num-1 );

                                      num_cofm_list++;
                                      cofm_list[num_cofm_list].x = c_of_m[0];
                                      cofm_list[num_cofm_list].y = c_of_m[1];
                                      cofm_list[num_cofm_list].z = c_of_m[2];
                                   }

                                 }
	/*********************************************/

                                  printf("\n");
                                  for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
				    {
                                      pdb_frame[num_pdb_atoms] = molecule[iloop];
 
                                      printf("atom: %s %10.6f %10.6f %10.6f\n",pdb_frame[num_pdb_atoms].label,
                                                                               pdb_frame[num_pdb_atoms].x,
                                                                               pdb_frame[num_pdb_atoms].y,
                                                                               pdb_frame[num_pdb_atoms].z);
 
                                      printf("image_disp %10.6f %10.6f %10.6f\n",image_disp[iloop].x,
                                                                                 image_disp[iloop].y,
                                                                                 image_disp[iloop].z);
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
			    printf("Molecule starts at %d ends at %d\n", start_mol, demarcation[index].end);

	/*********************************************/
	/* HAVE FLAG TO GIVE CO_ORDS *****************/
	/* DO NOT GO TO COFM IF NO MASSES SET  *******/
	/*********************************************/

				 printf("Atoms:\n");
				 write_car( &header_line, &title_line, &date_line,
				 &molecule[start_mol], pbc, &abc[0], demarcation[index].num-1, FALSE);

				 for (iloop=start_mol; iloop <= demarcation[index].end; iloop++)
				    {
				       printf("%3d %8s %10.6f %10.6f %10.6f\n", iloop, molecule[iloop].label,
					     molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);
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
 
                            printf("\n\nTotal velocity of molecule = %10.6f, %10.6f, %10.6f.\n", vel_list[vel_icorr].x,
                                                                                                 vel_list[vel_icorr].y,
                                                                                                 vel_list[vel_icorr].z);

                            vel_icorr++;
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

                  write_csv(msd_fp, "time", "msd", "msd", 
                                 &time_dots[0] , &msd[0], &msd[0], FALSE, 
                                max_msd_ind_diff-min_msd_ind_diff);

                  fclose(correlation_fp);
               }

	     if (anal_flags[ianal].big_med_small[0]) 
	       {
		  dot_correlation(&la_moi_vs_time[0], &la_moi_dots[0],
                          &time_dots[0], time_step, min_corr_ind_diff, 
                          max_corr_ind_diff,
                          frame_index, num_this_mol[1]-1 ); 

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
                  dot_correlation(&vel_list[0], &vel_dots[0],
                                  &time_dots[0], time_step,
                                  min_corr_ind_diff, max_corr_ind_diff,
                                  frame_index, num_this_mol[1]-1 );

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
