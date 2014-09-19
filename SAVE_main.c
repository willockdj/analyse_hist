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

double size_vector(double *p_vector);

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol );

void moments_of_inertia(atom *p_molecule, int num_atoms, double *p_c_of_m,
                           double *p_m_of_inertia, double *p_eigenvals );

void radius_gyration(atom *p_molecule, int num_atoms, double *p_c_of_m,   
                     double total_mass, double *p_rgyr );

void eigen_vec_3b3( double *p_matrix, double eigen_val, double *p_eigen_vec );

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, 
                double *p_recip_latt_vec, coord_flags *p_fix_flags,
                int need_draw, win_details *p_win_geom, int num_windows );

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
                types *p_monit_bnds, int *p_num_std_bnds, 
                monitors *p_monitor_set,  
                monitors *p_monitor_set_angle, 
                monitors *p_monitor_set_dihedral, 
                int *p_num_monit,
                int *p_num_monit_angle,
                int *p_num_monit_dihedral,
                int *p_num_limit,
                int *p_num_limit_angle,
                int *p_num_limit_dihedral,
                char *p_xdatcar_file, int *p_have_history, int *p_have_xdatcar, char *p_poscar_file, 
                char *p_potcar_file, int *p_have_poscar, int *p_have_potcar, windef *p_winref, 
                int *p_need_draw, double *p_rgyr_av_start);

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

void write_distrib_csv(FILE *fp, char *p_title_x, char *p_title_y, 
                       char *p_title_z,
                       double *p_x, int *p_y, int *p_z,
                       int have_z, int num);

void write_window_csv(FILE *fp, double *p_time_now, win_geometry *p_win_rads,
                      int start_type, int num_mols, int num_wins, int last_in_frame);

void write_window_csv_titles(FILE *fp,
                             int num_mols, int num_wins);

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

int find_window(atom *p_molecule, int num_atoms, windef *p_winref, winlist *p_winsets );

double circumcircle(atom *p_molecule, winlist *p_winsets, int iwin, double *p_centre, 
                    double *p_norm, double *p_rvec);

int main (int argc, char **argv)
{
FILE *fptr;
FILE *history_fp;
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
char xdatcar_file[80];
char poscar_file[80];
char potcar_file[80];
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
int written_coords=FALSE;
int written_windows=FALSE;

FILE *field_fp;
FILE *rdf_fp;
FILE *pdb_fp;
FILE *car_fp;
FILE *out_fp;
FILE *stat_fp;
FILE *msd_fp;
FILE *windows_fp;
FILE *err_fp;
FILE *correlation_fp;

int have_field, levcfg, traj_num;
int nstep, ianal, hist_state; 

double tstep[MAXFRAMES];
double time_now[MAXFRAMES];
double msd[MAXFRAMES];
double time_dots[MAXFRAMES];
double c_of_m[3], total_mass;
double m_of_inertia[6], moi_eigenvals[3], moi_eigenvecs[9];
double angle[3];
double big_eigen, scale, vec[3], vec1[3], vec2[3], vec3[3], temp, abc[6];
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
int num_bins, ineigh, this_neigh, neigh_index;
int iii, ibin, iguess, bins[MAX_BINS], bins_angle[MAX_BINS];
int iwin, iwin_tot, next_win, found_win, index_win;

/***********************************************/
/**** Updates for arc file writing, Sept 09 ****/
/***********************************************/
coord_flags fix_flags[MAX_ATOMS]; 
int start_frame, start_type;

/***********************************************/
/** Variables for H-bond monitoring ************/
/** Chris Lee Nov. 2010             ************/
/***********************************************/

              #define MAX_FILENAME 20
              int iother, bin_index;
              int iatom1, iatom2, iatom3, iatom4;
              int is_atom1, is_atom2, is_atom3, is_atom4, do_other;
              int have_history, have_xdatcar, have_poscar, have_potcar;
              double x1, x2, y1, y2, z1, z2 = 0;
              double x_atom1, x_atom2, x_atom3, x_atom4 = 0;
              double y_atom1, y_atom2, y_atom3, y_atom4 = 0;
              double z_atom1, z_atom2, z_atom3, z_atom4 = 0;
              double xa, xb, ya, yb, za, zb, dot, theta;
              int    num_angles, num_atom1;
              int    num_dihedrals = 0;
              double distance = 0;
              double distance_1_2 = 0; 
              double distance_2_3 = 0;
              double distance_3_4 = 0; 
              double bond_angle[MAX_MONIT];
              char fname_out[MAX_FILENAME];
              double four_pi_delr[MAX_MONIT];
              int num_sampled[MAX_MONIT]; /*** for number of pairs tested ***/
              int num_accepted[MAX_MONIT]; /*** for number of pairs accepted in distribution ***/
              int num_accepted_angles[MAX_MONIT]; /*** for number of pairs accepted in distribution ***/
              int num_accepted_dihedrals[MAX_MONIT]; /*** for number of pairs accepted in distribution ***/
              double cumulat_dist[MAX_MONIT];
              double bin_width[MAX_MONIT];
              double bin_width_angle[MAX_MONIT];
              double bin_width_dihedral[MAX_MONIT];
              monitors monitor_set[MAX_MONIT];
              monitors monitor_set_angle[MAX_MONIT];
              monitors monitor_set_dihedral[MAX_MONIT];
              double distrib[MAX_MONIT][MAX_BINS];
              double distrib_angle[MAX_MONIT][MAX_BINS];
              double distrib_dihedral[MAX_MONIT][MAX_BINS];
              double bin_distances[MAX_MONIT][MAX_BINS];
              double bin_angle[MAX_MONIT][MAX_BINS];
              double bin_dihedral[MAX_MONIT][MAX_BINS];
              double mean_dist[MAX_MONIT];
              int num_monit;
              int num_monit_angle;
              int num_monit_dihedral;
              int num_limit_bond;
              int num_limit_angle;
              int num_limit_dihedral;
              double vec_ux, vec_uy, vec_uz, vec_vx, vec_vy, vec_vz;
              double mod_u, mod_v;
              double dot_uv = 0;
              double phi = 0;

/*****************************************************************/
/** Variables for radius of gyration and for Windows analysis ****/
/** Dave Willock, March 2011.                                 ****/
/**                                                           ****/
/** rgyr_molav will contain the average rgyr for all molecules****/
/** in the indexed frame.                                     ****/
/*****************************************************************/

windef winref;
winlist winsets;
win_geometry win_rads;
win_details  win_geom;
double rgyr_molav[MAXFRAMES];
double mol_rgyr[MAXMOL];
double rgyr_av_start;
double cir_centre[3], cir_norm[3], cir_rvec[3], rad;
int need_draw, first_of_frame, last_in_frame;
int first_call;

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
printf("Going to read input\n");

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
                       &monit_bnds[0], &num_std_bnds, 
                       &monitor_set[0],  
		       &monitor_set_angle[0], 
                       &monitor_set_dihedral[0], 
		       &num_monit, 
		       &num_monit_angle, 
		       &num_monit_dihedral, 
                       &num_limit_bond,
                       &num_limit_angle, 
                       &num_limit_dihedral,
                       &xdatcar_file[0], &have_history, &have_xdatcar,
                       &poscar_file[0], &potcar_file[0], &have_poscar, &have_potcar,
                       &winref, &need_draw, &rgyr_av_start );

printf("Going to read input back\n");
/**** Check file type ***/

if (have_xdatcar && have_history)
  {
    printf("ERROR: Both HISTORY and XDATCAR files defined, only give one trajectory file!\n");
    exit(0);
  }
else if  (!have_xdatcar && !have_history && !have_config)
  {
    printf("ERROR: None of HISTORY, XDATCAR or CONFIG files are defined, need a trajectory or structure file!\n");
    exit(0);
  }
else if (have_xdatcar)
  {
    printf("This is an XDATCAR file, processing VASP trajectory\n");
    if (have_potcar)
      {
         printf("Also have POTCAR file defined\n");
      }
    else
      {
         printf("ERROR: XDATCAR file supplied without corresponding POTCAR file!\n");
         exit(0);
      }
  }
else if (have_history)
  {
    printf("This is a HISTORY file, processing dlpoly trajectory\n");
  }
else if (have_config)
  {
    printf("This is a CONFIG file, processing dlpoly input structure.\n");
  }
else if (have_poscar)
  {
    printf("This is a POSCAR file, prosessing VASP input structure.\n");
  }

/*** Note num_monit is the number of monitors set and num_limits the number of ranges defined ***/
/*** Need to add in a test for these are the same if monitors required.                       ***/

printf("Have num_monit=%d\n", num_monit);
if (num_monit >= 0)
  {
     if  (num_limit_bond != num_monit)
       {
	  printf("ERROR: The number of bonds set to monitor does not equal the number of limits defined\n");
	  exit(0);
       }

     printf("Requested monitoring of atom-atom interactions:\n");
     for (iii=0; iii<=num_monit; iii++)
       {
         printf("%4s %4s in range %10.6f - %10.6f with %d bins\n", monitor_set[iii].atom1,
                                                                   monitor_set[iii].atom2,
                                                                   monitor_set[iii].min_limit,
                                                                 monitor_set[iii].max_limit,
                                                                 monitor_set[iii].bins);  
                                                                 
         if (monitor_set[iii].bins >= MAX_BINS)
           {
              printf("ERROR: That monitor set has too many bins defined current maximum is : %d\n", MAX_BINS); 
              exit(0);
           }
       }
  }

printf("Have num_monit_angle=%d\n", num_monit_angle);
if (num_monit_angle >= 0)
  {
     if  (num_limit_angle != num_monit_angle)
       {
	  printf("ERROR: The number of angles set to monitor does not equal the number of limits defined\n");
	  exit(0);
       }
     printf("Requested monitoring of bond angles:\n");
     for (iii=0; iii<=num_monit_angle; iii++)
       {
         printf("%4s %4s %4s with %d bins\n", monitor_set_angle[iii].atom1,
                                              monitor_set_angle[iii].atom2,
                                              monitor_set_angle[iii].atom3,
                                              monitor_set_angle[iii].bins);

         printf("Minimum angle to record = %10.6f\n", monitor_set_angle[iii].min_limit);
         printf("Maximum angle to record = %10.6f\n", monitor_set_angle[iii].max_limit);

         printf("Limit for 1-2 = %10.6f\n", monitor_set_angle[iii].max_limit1);
         printf("Limit for 2-3 = %10.6f\n", monitor_set_angle[iii].max_limit2);

         if (monitor_set_angle[iii].bins >= MAX_BINS)
           {
              printf("ERROR: That monitor set has too many bins defined current maximum is : %d\n", MAX_BINS);
              exit(0);
           }
       }
  }

printf("Have num_monit_dihedral=%d\n", num_monit_dihedral);
if (num_monit_dihedral >= 0)
  {
     if  (num_limit_dihedral != num_monit_dihedral)
       {
	  printf("ERROR: The number of dihedral angles set to monitor does not equal the number of limits defined\n");
	  exit(0);
       }
     printf("Requested monitoring of dihedral angles:\n");
     for (iii=0; iii<=num_monit_dihedral; iii++)
       {
         printf("%4s %4s %4s %4s with %d bins\n", monitor_set_dihedral[iii].atom1,
                                                  monitor_set_dihedral[iii].atom2,
                                                  monitor_set_dihedral[iii].atom3,
                                                  monitor_set_dihedral[iii].atom4,
                                                  monitor_set_dihedral[iii].bins);

         printf("Minimum angle to record = %10.6f\n", monitor_set_dihedral[iii].min_limit);
         printf("Maximum angle to record = %10.6f\n", monitor_set_dihedral[iii].max_limit);

         printf("Limit for 1-2 = %10.6f\n", monitor_set_dihedral[iii].max_limit1);
         printf("Limit for 2-3 = %10.6f\n", monitor_set_dihedral[iii].max_limit2);
         printf("Limit for 3-4 = %10.6f\n", monitor_set_dihedral[iii].max_limit3);

         if (monitor_set_dihedral[iii].bins >= MAX_BINS)
           {
              printf("ERROR: That monitor set has too many bins defined current maximum is : %d\n",
                                                                                           MAX_BINS);
              exit(0);
           }
       }
  }

printf("Back from read_input with out_file: >>%s<<\n", out_file);
if (DEBUG) printf("DEBUG mode is on\n");
else  printf("DEBUG mode is off\n");

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
     
     printf("FIELD file says there are %d different types of molecule\n", num_types);
     for (iloop=0; iloop <= num_types; iloop++)
        {
           printf("Molecule %d has %d instances\n", iloop, num_this_mol[iloop]);
           printf("....starts %d and has %d atoms\n", demarcation[iloop].start,  demarcation[iloop].num);
        }
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
    else
     {
        num_stat_list=-1;
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

printf("DLPOLY run used a time between frames of %10.6f ps\n", time_step );

printf("%d Analysis passes requested\n", num_anal_flags+1);

for (ianal=0; ianal <= num_anal_flags; ianal++)
  {
    printf("\nAnalysis pass %d is for:\n", ianal+1);
    if (anal_flags[ianal].mo_inertia) printf("Moment of inertia calculation\n");
    if (anal_flags[ianal].big_med_small[0])
      {
         printf("Largest eigenvalue will be reported\n");
         if (anal_flags[ianal].sum_comps)
              printf("Totals and averages of the unit vector components will also be generated\n");
      }
    if (anal_flags[ianal].big_med_small[1])
      {
         printf("Middle eigenvalue will be reported\n");
         if (anal_flags[ianal].sum_comps)
              printf("Totals and averages of the unit vector components will also be generated\n");
      }
    if (anal_flags[ianal].big_med_small[2])
      {
         printf("Smallest eigenvalue will be reported\n");
         if (anal_flags[ianal].sum_comps)
              printf("Totals and averages of the unit vector components will also be generated\n");
      }
    if (anal_flags[ianal].correlation)
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
    if (anal_flags[ianal].msd)
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
    if (anal_flags[ianal].angles)
                   printf("Will report angles between MOI vector and axes\n");
    if (anal_flags[ianal].dipoles)
                   printf("Will report molecular dipole moments\n");
    if (anal_flags[ianal].forster)
                   printf("Will report forster dipole orientation factor averaged over each frame\n");
    if (anal_flags[ianal].coords)
       {
          printf("Will print molecular co-ordinates for each molecule so that all atoms are\n");
          printf("in a common unit cell.\n");
       }
    if (anal_flags[ianal].contact)
      {
        printf("Will search for all inter-molecular close contacts below %10.6f Angstroms\n",
                cutoff);
        cutoff2= cutoff*cutoff;
      }
    if (anal_flags[ianal].add_atom)
      {
        printf("Will add atoms at the geometric positions defined in the input file\n");
        if (anal_flags[ianal].bisector )
          {
             printf("A new site will be added along the bisector of %s %s %s triplets of atoms", 
                           bisector_atoms[0].label, bisector_atoms[1].label, bisector_atoms[2].label);
             printf(" %10.6f from the central atom into the acute angle\n", bisect_disp);
          }
      }
    if (anal_flags[ianal].pdb)
      {
        printf("Will output trajectory file in pdb format as %s\n", pdb_file);
        if ( want_raw )
          {
             printf("Will print raw pdb data, i.e. with periodic boundaries shown\n");
          }
      }
   
   if (anal_flags[ianal].bonds)
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

/*** If required open file for windows ****/
   if (anal_flags[ianal].windows)
     {
        sprintf(fname_out, "windows.csv");

        if ( (windows_fp = fopen (&fname_out[0],"w") ) == NULL)
          {
             printf ("ERROR Opening %s FILE\n", fname_out);
             return 1;
          }
        else written_windows=TRUE;
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

	sm_icorr = 0;
	me_icorr = 0;
	la_icorr = 0;
	abc_form= FALSE;
	frame_index = -1;
        num_cofm_list= -1;

/*************************************************************************/
/*** Set up bins for H-bond monitoring ***********************************/
/*************************************************************************/
     if (num_monit >= 0)
       {
        for (iii=0; iii<=num_monit; iii++)
           {
                num_sampled[iii]=0;
                num_accepted[iii]=0;
                bin_width[iii]= ( monitor_set[iii].max_limit - monitor_set[iii].min_limit)/monitor_set[iii].bins;
                
/*** First part of normalisation factor for g(r) plots *******************/

                four_pi_delr[iii]= four_pi * bin_width[iii];

                for (icheck=0; icheck < monitor_set[iii].bins; icheck++)
                   {
                      distrib[iii][icheck]=0;
                      bin_distances[iii][icheck]=monitor_set[iii].min_limit+(icheck+0.5)*bin_width[iii];
                   }
           }
       }
     if (num_monit_angle >= 0)
       {
        for (iii=0; iii<=num_monit_angle; iii++)
           {
                num_accepted_angles[iii]=0;
                bin_width_angle[iii] = ( monitor_set_angle[iii].max_limit 
                                       - monitor_set_angle[iii].min_limit)
                                                      /monitor_set_angle[iii].bins;
                
                for (icheck=0; icheck < monitor_set_angle[iii].bins; icheck++)
                   {
                    distrib_angle[iii][icheck]=0;
                    bin_angle[iii][icheck]= monitor_set_angle[iii].min_limit +
                                                       (icheck+0.5)*bin_width_angle[iii];

                    if (DEBUG) printf("Set bin_angle %d for case %d to %10.6f\n", 
                                                   icheck, iii, bin_angle[iii][icheck]);
                   }
           }
       }
     if (num_monit_dihedral >= 0)
       {
        for (iii=0; iii<=num_monit_dihedral; iii++)
           {
                num_accepted_dihedrals[iii]=0;
                bin_width_dihedral[iii] = ( monitor_set_dihedral[iii].max_limit 
                                          - monitor_set_dihedral[iii].min_limit)
                                                      /monitor_set_dihedral[iii].bins;
                
                for (icheck=0; icheck < monitor_set_dihedral[iii].bins; icheck++)
                   {
                    distrib_dihedral[iii][icheck]=0; 
                    bin_dihedral[iii][icheck]= monitor_set_dihedral[iii].min_limit 
                                                       +(icheck+0.5)*bin_width_dihedral[iii];

                    if (DEBUG) printf("Set bin_dihedral %d for case %d to %10.6f\n", 
                                                 icheck, iii, bin_dihedral[iii][icheck]);
                   }
           }
       }
        num_angles=0;
        num_atom1=0;

	/*********************************************************************************/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/***** Read a frame and process it, frame reader returns -10 at end of file ! ****/
	/*********************************************************************************/

	while ( read_hist_frame( history_fp, &molecule[0], max_atom_index, &nstep, 
                                 &tstep[frame_index+1], 
			         pbc, &latt_vec[0], have_config, levcfg ) > 0) 
	  {
	      frame_index++;

              if (frame_index >= MAXFRAMES)
                {
                   printf("ERROR : Too many frames in HISTORY file for current version.\n");
                   printf("        Recompile with MAXFRAMES increased.\n");
                   exit(0);
                }
              if (frame_index == 0) start_frame=TRUE;
                               else start_frame=FALSE;

              time_now[frame_index]= nstep*tstep[frame_index];
              printf("Processing frame %d timestep: %10.6f time: %10.6f ps\n",frame_index, 
                                                                              time_now[frame_index],
                                                                              tstep[frame_index]);
              first_of_frame=TRUE;

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

/**** Use information from the FIELD file to do neighbours on a molecule by molecule basis. ***/

              first_call=TRUE;
              for (iloop=0; iloop <= num_types; iloop++)
                 {
                     start_mol= demarcation[iloop].start;

                     generate_neighbours( &molecule[start_mol], demarcation[iloop].num-1, 
                                          &atom_types[0], &num_atom_types,
                                          TRUE, first_call);
                     first_call = FALSE;
                 }

              if (start_frame || DEBUG)
                {
                  printf("\nNeighbour information for frame %d:\n", frame_index);
                  index=0;
                  for (iloop=0; iloop <= num_types; iloop++)
                    {
                     start_mol= demarcation[iloop].start;
                     printf("\nMolecule %d with %d atoms.......\n\n",iloop, demarcation[iloop].num);

                     for (icheck=0; icheck < demarcation[iloop].num; icheck++)
                       {
                         printf("%d (abs: %d), %s (%s) : ", icheck, index,
                                                    molecule[index].label,
                                                    molecule[index].elem);

                         for (ineigh=0; ineigh<= molecule[index].num_neigh; ineigh++)
                            {
                               neigh_index=start_mol+molecule[index].neighb[ineigh];

                               vec[0]= molecule[index].x - molecule[neigh_index].x;
                               vec[1]= molecule[index].y - molecule[neigh_index].y;
                               vec[2]= molecule[index].z - molecule[neigh_index].z;

                               min_image(&vec[0], &vec[1], &vec[2]);

                               printf(" %s (%s, d=%6.3f)", molecule[neigh_index].label,
                                                           molecule[neigh_index].elem,
                                                           size_vector(&vec[0]));
                            }
                         printf("\n");
                         index++;
                      }
                    }
                }

/** Start of H-bond processing additions                              */
/*put all labels and corresponding distances into arrays              */
/*calculate and print average distances for h-bonding for this frame  */
/*move to next frame and repeat calculation                           */
/*get averages across all frames and calculate global average         */
/*Chris Lee Dec. 2010                                                 */

/*************************************************************************/
/*** Process frame for atom-atom distance g(r) monitoring  ***************/
/***                   angle and dihedral histograms.      ***************/
/*** These require no molecule information and so are done ***************/
/*** without regard to molecule indexes or to the use_this ***************/
/*** array setting.                                        ***************/
/*** Dave Willock July 2011                                ***************/
/*************************************************************************/

/*************************************************************************/
/*** Set up bins for H-bond monitoring ***********************************/
/*************************************************************************/
      for (iii=0; iii<=num_monit; iii++)
         {

          for (icheck=0; icheck<=max_atom_index; icheck++)                  
              {
                  is_atom1 = FALSE;
                  is_atom2 = FALSE;

                if (strcmp(molecule[icheck].label,monitor_set[iii].atom1) == 0)
                  {
                    is_atom1 = TRUE;
                  }
                else if (strcmp(molecule[icheck].label,monitor_set[iii].atom2) == 0)
                  {
                    is_atom2 = TRUE;
                  }

                if (is_atom1 || is_atom2)
                  {
                    if (DEBUG)
                      {
                        printf("Checking label %d >>%s<< elem >>%s<< x: %10.6f y: %10.6f z: %10.6f\n", icheck,               
                                                          molecule[icheck].label,
                                                          molecule[icheck].elem,
                                                          molecule[icheck].x,
                                                          molecule[icheck].y,
                                                          molecule[icheck].z );
                      }
             /*begin hydrogen bonding calculation*/             

                     x1=  molecule[icheck].x;
                     y1=  molecule[icheck].y;
                     z1=  molecule[icheck].z;
             
                     for (iother=icheck+1; iother<=max_atom_index; iother++)                  
                        {

             /*start with first frame*/ 
             /* for each pair of atoms with appropriate labels */
                
                if (( is_atom1 && strcmp(molecule[iother].label,monitor_set[iii].atom2) == 0)
                    || (strcmp(molecule[iother].label,monitor_set[iii].atom1) == 0 && is_atom2))
                  {
                    num_sampled[iii]++;

                    x2=  molecule[iother].x;
                    y2=  molecule[iother].y;
                    z2=  molecule[iother].z;

                    vec[0] = x1-x2;
                    vec[1] = y1-y2;
                    vec[2] = z1-z2;

                    min_image(&vec[0], &vec[1], &vec[2]);

                    distance = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                     
                     if (monitor_set[iii].min_limit <= distance && distance < monitor_set[iii].max_limit)            
                       {
                        if (DEBUG)
                           {
                             printf("DEBUG>> Accepted pair of atoms in range %10.6f to %10.6f in bins with width %10.6f\n",
                                                               monitor_set[iii].min_limit, monitor_set[iii].max_limit, 
                                                               bin_width[iii]);
                             printf("%s %s %10.6f\n", molecule[icheck].elem, molecule[iother].elem, distance);
                           }
                        num_accepted[iii]++;
                        cumulat_dist[iii] = cumulat_dist[iii] + distance;

                        bin_index = (distance-monitor_set[iii].min_limit)/bin_width[iii];
                        distrib[iii][bin_index]++;
                       }
                  }
                 }
               }
              }
        printf("Cumulative %s..%s distance = %10.6f\n", monitor_set[iii].atom1, 
                                                        monitor_set[iii].atom2,  
                                                        cumulat_dist[iii]); 
        }
     /*end hydrogen bonding calculation*/

   /****begin angle calculation for 1-2-3 ****/
   for (iii=0; iii<=num_monit_angle; iii++) 
      {
       printf("Calculating specified bond angles\n");
/*** Look for central atom ***/
       for (iatom2=0; iatom2<=max_atom_index; iatom2++)
          {
            if (strcmp(molecule[iatom2].label,monitor_set_angle[iii].atom2) == 0)
                  {
                     x_atom2=  molecule[iatom2].x;
                     y_atom2=  molecule[iatom2].y;
                     z_atom2=  molecule[iatom2].z;
 
                for (iatom1=0; iatom1<=max_atom_index; iatom1++)
                   {
                    if (strcmp(molecule[iatom1].label,monitor_set_angle[iii].atom1) == 0)   
                          {
                           x_atom1=  molecule[iatom1].x;
                           y_atom1=  molecule[iatom1].y;
                           z_atom1=  molecule[iatom1].z;
                            
                           if (iatom2 != iatom1)
                             {
/*** test if iatom1 and iatom2 are close enough according to the limits set ***/
/*** No need for distance_1_2 to be an array                                ***/

                            vec[0] = x_atom1-x_atom2;
                            vec[1] = y_atom1-y_atom2;
                            vec[2] = z_atom1-z_atom2;

                            min_image(&vec[0], &vec[1], &vec[2]);

                            distance_1_2 = sqrt( vec[0]*vec[0]+ vec[1]*vec[1]+ vec[2]*vec[2]);

                            if (distance_1_2 < monitor_set_angle[iii].max_limit1)
                              {
                                num_atom1++;
/**** Look for third atom to make up angle ***/

                                for (iatom3=iatom1+1; iatom3<=max_atom_index; iatom3++)
                                   {
                                    if (strcmp(molecule[iatom3].label,monitor_set_angle[iii].atom3) == 0)
                                      {
                                        if (iatom3 != iatom2)
                                          {
                                            x_atom3=  molecule[iatom3].x;
                                            y_atom3=  molecule[iatom3].y;
                                            z_atom3=  molecule[iatom3].z;

                                            vec1[0] = x_atom3-x_atom2;
                                            vec1[1] = y_atom3-y_atom2;
                                            vec1[2] = z_atom3-z_atom2;

                                            min_image(&vec1[0], &vec1[1], &vec1[2]);

                                            distance_2_3 = sqrt( vec1[0]*vec1[0]+ vec1[1]*vec1[1]+ vec1[2]*vec1[2]);

                                            if (distance_2_3 < monitor_set_angle[iii].max_limit2)
                                              {
/*** define bond vectors based on 2 being the central atom, i.e. 2->1 and 2->3 ***/
                                                 if (DEBUG)
                                                      printf("Getting angle for index set %d %d %d\n", iatom1, iatom2, iatom3);

                                                 dot = vec[0]*vec1[0]+vec[1]*vec1[1]+vec[2]*vec1[2];

                                                 theta = RAD_TO_DEG * acos( dot / (distance_2_3 * distance_1_2));
/*** Increment the corresponding bin ***/
                                                 bin_index = (theta-monitor_set_angle[iii].min_limit)
                                                                                     /bin_width_angle[iii];

                                                 distrib_angle[iii][bin_index] += 1.0;

                                                 if (DEBUG)
                                                 printf("Angle mon %d: theta %10.6f -> bin %d centred %10.6f width %10.6f\n",
                                                                  iii, theta, bin_index, 
                                                                  bin_angle[iii][bin_index], bin_width_angle[iii] );  

                                                 num_accepted_angles[iii]++; 

                                              }
                                          }
                                       }
                                     }
                              }
                          }
                       }
                    }
/*** report failure ***/
                  }
                 }
               } 
               if (num_monit >= 0) printf("Number of angles found %d\n", num_angles);
   /****end angle calculation****/
  
 /****begin diheral angle calculation***/
  for (iii=0; iii<=num_monit_dihedral; iii++)
      {
        if (DEBUG)
        printf("Calculating specified dihedral angle %s %s %s %s\n", monitor_set_dihedral[iii].atom1
                                                                   , monitor_set_dihedral[iii].atom2
                                                                   , monitor_set_dihedral[iii].atom3
                                                                   , monitor_set_dihedral[iii].atom4);
        for (iatom1=0; iatom1<=max_atom_index; iatom1++)
          {

/** dihedral defined as 1-2-3-4 ***/
/** So looking for planes 1-2-3 ***/
/**                   and 2-3-4 ***/

           is_atom1 = FALSE;
           is_atom2 = FALSE;
           is_atom3 = FALSE;
           is_atom4 = FALSE;

           if (strcmp(molecule[iatom1].label,monitor_set_dihedral[iii].atom1) == 0)
             {              
                x_atom1 = molecule[iatom1].x;
                y_atom1 = molecule[iatom1].y;
                z_atom1 = molecule[iatom1].z;

                for (iatom2=0; iatom2<=max_atom_index; iatom2++)
                   {
                    if (strcmp(molecule[iatom2].label,monitor_set_dihedral[iii].atom2) == 0)
                      {
                       x_atom2 = molecule[iatom2].x;
                       y_atom2 = molecule[iatom2].y;
                       z_atom2 = molecule[iatom2].z;

                       if (iatom1 != iatom2)
                         {
                          vec[0] = x_atom1-x_atom2;
                          vec[1] = y_atom1-y_atom2;
                          vec[2] = z_atom1-z_atom2;

                          min_image(&vec[0], &vec[1], &vec[2]);

                          distance_1_2 = sqrt( vec[0]*vec[0]+ vec[1]*vec[1]+ vec[2]*vec[2]);

                          if (distance_1_2 < monitor_set_dihedral[iii].max_limit1)
                            {
                             for (iatom3=0; iatom3<=max_atom_index; iatom3++)
                                {
                                 if (strcmp(molecule[iatom3].label,monitor_set_dihedral[iii].atom3) == 0)
                                   {
                                    x_atom3 = molecule[iatom3].x;
                                    y_atom3 = molecule[iatom3].y;
                                    z_atom3 = molecule[iatom3].z;

                                    if (iatom3 != iatom2 && iatom3 != iatom1)
                                      {
                                       vec1[0] = x_atom2-x_atom3;
                                       vec1[1] = y_atom2-y_atom3;
                                       vec1[2] = z_atom2-z_atom3;
                                       
                                       min_image(&vec1[0], &vec1[1], &vec1[2]);

                                       distance_2_3 = sqrt( vec1[0]*vec1[0]+ vec1[1]*vec1[1]+ vec1[2]*vec1[2]);

                                       if (distance_2_3 < monitor_set_dihedral[iii].max_limit2)
                                         {
                                          for (iatom4=iatom1+1; iatom4<=max_atom_index; iatom4++)
                                             {
                                              if (strcmp(molecule[iatom4].label,monitor_set_dihedral[iii].atom4) == 0)
                                                {
                                                 x_atom4 = molecule[iatom4].x;
                                                 y_atom4 = molecule[iatom4].y;
                                                 z_atom4 = molecule[iatom4].z;

                                                 if (iatom4 != iatom3 && iatom4 != iatom2 && iatom4 != iatom1)
                                                   {
                                                    /*count number of dihedral angles found*/
                                                    num_dihedrals++;
                                                     
                                                    vec2[0] = x_atom3-x_atom4;
                                                    vec2[1] = y_atom3-y_atom4;
                                                    vec2[2] = z_atom3-z_atom4;

                                                    min_image(&vec2[0], &vec2[1], &vec2[2]);

                                                    distance_3_4 = sqrt( vec2[0]*vec2[0]+ vec2[1]*vec2[1]+ vec2[2]*vec2[2]);
                                                    
                                                    if (distance_3_4 < monitor_set_dihedral[iii].max_limit3)
                                                      {
                                                      /*do dihedral angle calculation*/
                                                       vec_ux =   vec[1]*vec1[2] - vec1[1]*vec[2];
                                                       vec_uy = -(vec[0]*vec1[2] - vec1[0]*vec[2]);
                                                       vec_uz =   vec[0]*vec1[1] - vec1[0]*vec[1];

                                                       vec_vx =   vec1[1]*vec2[2] - vec2[1]*vec1[2];
                                                       vec_vy = -(vec1[0]*vec2[2] - vec2[0]*vec1[2]);
                                                       vec_vz =   vec1[0]*vec2[1] - vec2[0]*vec1[1];  
                                                        
                                                       mod_u = sqrt(vec_ux*vec_ux + vec_uy*vec_uy + vec_uz*vec_uz);
                                                       mod_v = sqrt(vec_vx*vec_vx + vec_vy*vec_vy + vec_vz*vec_vz);

                                                       dot_uv = vec_ux*vec_vx + vec_uy*vec_vy + vec_uz*vec_vz;

                                                       phi = RAD_TO_DEG*acos((dot_uv/(mod_u*mod_v)));

                                                      /*sign of phi from position of 34 vector over 123 plane */

                                                       dot = -(vec_ux * vec2[0] +vec_uy * vec2[1] +vec_uz * vec2[2]);

                                                       if (dot < 0.0) phi=-phi;

                                                       if (DEBUG)
                                                         {
                                                       printf("Getting dihedral angle for index set %d %d %d %d\n",
                                                               iatom1, iatom2, iatom3, iatom4);
                                                       printf("Distances: 1-2 = %10.6f, 2-3 = %10.6f, 3-4 = %10.6f\n",
                                                               distance_1_2, distance_2_3, distance_3_4);

                                                       printf("labels: %s %s %s %s\n", molecule[iatom1].label,
                                                                                       molecule[iatom2].label, 
                                                                                       molecule[iatom3].label, 
                                                                                       molecule[iatom4].label);

                                                       printf("Dihedral angle = %10.6f\n", phi);
                                                       printf("Number of dihedral angles found = %d\n", num_dihedrals);
                                                         }

/*** Increment the corresponding bin remenbering dihedrals run -180.0 to 180.0 ***/

						       num_accepted_dihedrals[iii]++;
                                                       bin_index = (phi-monitor_set_dihedral[iii].min_limit)
                                                                                     /bin_width_dihedral[iii];

                                                       distrib_dihedral[iii][bin_index] += 1.0;
                                                      }
                                                   }
                                                }
                                             }
                                         }
                                      }
                                   }
                                }       
                            }
                         }
                      }
                   }
               }
   
          }       
      }  
   /***end diheral angle calculataion****/

/**********************************************************************/
/** Analysis passes for processes that rely on molecule indices *******/
/**********************************************************************/
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

        /*********************************************/
        /** Loop over the molecule types identified **/
        /*********************************************/

                  start_type=TRUE;
		  for (itype=0; itype < num_molecules; itype++)
		    {

                      last_in_frame=TRUE;
                      for (jtype = itype+1; jtype < num_molecules; jtype++) 
                                          if (use_type[jtype]) last_in_frame=FALSE;
                      
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

                                if (DEBUG) printf("Top of this_molecule loop with %d index %d\n", this_molecule, index);

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
                                  if (start_frame && this_molecule == 0 && !done_pdb)
                                    {
	                               if ( anal_flags[ianal].pdb && (pdb_fp = fopen (&pdb_file[0],"w") ) == NULL)
                                        {
                                         printf ("ERROR Opening pdb file %s for writing\n", pdb_file);
                                         perror ("open");
                                         return 1;
	                                }
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
				    			                    demarcation[index].num-1, -1 );

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

                                       min_image(&vec[0],&vec[1],&vec[2]);
				       unit_vector(&vec[0]);

				       vec1[0] = molecule[bis_index3].x - molecule[bis_index2].x;
				       vec1[1] = molecule[bis_index3].y - molecule[bis_index2].y;
				       vec1[2] = molecule[bis_index3].z - molecule[bis_index2].z;

                                       min_image(&vec[0],&vec[1],&vec[2]);
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
							       demarcation[index].num-1, -1 );
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
/*** Added for cage structures ****/
                      if (anal_flags[ianal].windows)
                         {

/****Note ONLY FIRST FRAME IS USED TO DEFINE WINDOWS                          ***/
/****We are in a loop over the molecules controlled by this_molecule variable ***/
                            if (start_frame)
                              {
                                if (first_of_frame)
                                  {
                                    iwin_tot = -1;
                                    next_win = 0;
                                    winsets.num_windows = -1;
                                  }
                                printf("Going to find windows\n");
                                printf("Expecting %d per molecule\n", winref.num_windows);
                                printf("Using %d atoms to define a window.\n", winref.num_atoms);
                                for (iatom=0; iatom < winref.num_atoms; iatom++)
                                   {
                                  printf("Ref: %d Atom in window: %s neighbouring %s\n", iatom,
                                                                          &(winref.atoms[iatom][0]),
                                                                          &(winref.neigh[iatom][0]));
                                   }

                                found_win= find_window(&molecule[start_mol], demarcation[index].num-1, 
                                                       &winref, &winsets );

                                if (found_win)
                                  {
                                    printf("\nFound %d windows so far.\n", winsets.num_windows+1 );

                                    if (winsets.num_windows - next_win != winref.num_windows-1)
                                       {
                                          printf("ERROR: molecule %d does not contain the expected number of windows,\n", index);
                                          printf("ERROR: it actually has %d check use_type setting\n", winsets.num_windows - next_win +1 );
                                          printf("ERROR: Also check structure in error.car\n" );
                                          sprintf(fname_out, "error.car");

	                                  if ( (err_fp = fopen (&fname_out[0],"w" )) == NULL)
                                            {
                                               printf ("ANOTHER ERROR Opening error.car for writing!\n");
                                               exit(0);
	                                    }

                                         printf("Going to write_car with %s\n", fname_out);
                                         write_car( err_fp, &header_line[0], &title_line[0], &c_title_line[0],
		                                    &date_line[0], &molecule[start_mol], &mol_nums[0], 
                                                    pbc, &abc[0], demarcation[index].num, 1.0, 
                                                    start_frame, &super[0], &latt_vec[0], &recip_latt_vec[0], 
                                                    &fix_flags[0], need_draw, &win_geom, iwin_tot);

                                         exit(0);

                                       }
                                    for (iwin=next_win; iwin <= winsets.num_windows; iwin++)
                                      {
                                        iwin_tot++;
                                        if (iwin_tot >= MAX_WIN_TOT) 
                                          {
                                             printf("ERROR: number of windows now exceeds maximum allowed by MAX_WIN_TOT ( %d )\n", MAX_WIN_TOT);
                                             exit(0);
                                          }
                                        printf("Window %d defined by atoms; ", iwin);
                                        for (iii=0; iii<= winsets.num_atoms[iwin]; iii++)
                                          {
                                            iatom=winsets.iatom[iwin][iii];
                                            printf("%d (%s) ", iatom, molecule[start_mol+iatom].label);
                                            sprintf(molecule[start_mol+iatom].label,"%s%d","W",iwin);
                                          }
                                       printf("\n");
                                     }
                                    next_win=iwin_tot+1;
                                    printf("Copied over to next_win now %d\n", next_win);
                                  }
                               }
/*** Obtain circumcirle centre, radius and normal for all windows in the molecule set ****/

                           for (iwin=0; iwin < winref.num_windows; iwin++)
                             {
                              index_win= this_molecule*(winsets.num_windows+1)+iwin;
                              win_rads.radius[index_win]= circumcircle(&molecule[start_mol], 
                                                                       &winsets, 
                                                                       iwin, &cir_centre[0],
                                                                       &cir_norm[0],
                                                                       &cir_rvec[0]);

                              for (iii=0; iii<3; iii++) 
                                {
                                    win_geom.centre[index_win][iii] = cir_centre[iii];
                                    win_geom.norm[index_win][iii] = cir_norm[iii];
                                    win_geom.r_vec[index_win][iii] = cir_rvec[iii];
                                }
                             }
                           }
 
/*** Added for radius of gyration ***/
                      if (anal_flags[ianal].rgyr)
                         {
			   centre_of_mass(&c_of_m[0], &total_mass, &molecule[start_mol], 
							       demarcation[index].num-1, -1 );

                           radius_gyration( &molecule[start_mol], demarcation[index].num-1, &c_of_m[0],
                                            total_mass, &mol_rgyr[this_molecule]);

                         }

       if (anal_flags[ianal].coords)
           {
	     if (DEBUG) printf("In coords processing: Frame %d, molecule %d starts at %d ends at %d\n", 
                                        frame_index, index, start_mol, demarcation[index].end);

	/*********************************************/
	/* HAVE FLAG TO GIVE CO_ORDS *****************/
	/* DO NOT GO TO COFM IF NO MASSES SET  *******/
	/*********************************************/

        /*** Just output all molecules when we get to the first one for complete frames in arc file ***/
            if (start_mol == 0)
              {
                     if (start_frame) 
                        {
                           printf("Openning movie.arc\n");
                           sprintf(car_file,"movie.arc");
                           sprintf(c_title_line,"frame number %d from dlpoly run\n",frame_index+1);

        /** flag header, date, and title as blank so defaults are used ***/
                           header_line[0]=-1;
                           title_line[0]=-1;
                           date_line[0]=-1;

        /** Assume 1 by 1 by 1 cell should be output ****/
                           for (iii=0; iii<3; iii++) super[iii] = 1;

                           if ( (car_fp = fopen (&car_file[0],"w") ) == NULL)
                             {
                                printf("ERROR openning car file for writing\n");
                                exit(0);
                             }
                           written_coords=TRUE;
                         }

                      for (iloop=0; iloop < max_atom_index; iloop++)
                         {
                            mol_nums[iloop] = 0;
                            fix_flags[iloop].fx=FALSE;
                            fix_flags[iloop].fy=FALSE;
                            fix_flags[iloop].fz=FALSE;
                            molecule[iloop].mol=0;
                         }

/*                               printf("...gathering..with max_atom_index=%d.\n",max_atom_index); */
/* For gather molecule should do one molecule at a time not the entire set ****/
/*                                 gather_molecule(&molecule[start_mol], max_atom_index, 0); */

                       write_car( car_fp, &header_line[0], &title_line[0], &c_title_line[0],
		                  &date_line[0], &molecule[start_mol], &mol_nums[0], 
                                  pbc, &abc[0], max_atom_index+1, 1.0, 
                                  start_frame, &super[0], &latt_vec[0], &recip_latt_vec[0], 
                                  &fix_flags[0], need_draw, &win_geom, iwin_tot);

                         if (DEBUG)
                           {
                             for (iloop=start_mol; iloop <= max_atom_index; iloop++)
			      {
			       printf("%3d %8s %10.6f %10.6f %10.6f\n", iloop, molecule[iloop].label,
				     molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);
			      }
                           }
                         }
                 if (DEBUG) printf("....Done coords processing...\n");
	       }

	      index++;

/**** end of this_molecule loop but NOT itype ***/
           }

/**** Write window information for this set of molecules now ****/
             if (anal_flags[ianal].windows)
               {

                  if (start_frame && start_type)
                    {
                         write_window_csv_titles(windows_fp,
                                                 num_molecules, winref.num_windows);
                    }
        /*********************************************/
        /*** Only include flagged molecule types   ***/
        /*********************************************/
                  for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
		    {
                       printf("Writing windows for molecule type %d has %d instances and %d windows\n",
                                itype, num_this_mol[itype], winref.num_windows);

                       write_window_csv(windows_fp, &time_now[frame_index], &win_rads,
                                        start_type, num_this_mol[itype], 
                                        winref.num_windows, last_in_frame);
                   }
               }

/*** Average values across molecules for this frame as required ****/
                   if (anal_flags[ianal].rgyr) rgyr_molav[frame_index] = 0.0;

		   for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
   	              {
                        if (anal_flags[ianal].rgyr) 
                                   rgyr_molav[frame_index] += mol_rgyr[this_molecule];
		      }

                   if (anal_flags[ianal].rgyr)
                     {
                         rgyr_molav[frame_index] = rgyr_molav[frame_index]/num_this_mol[itype];

                         printf("Average radius of gyration across molecules for frame %d = %10.6f Angstroms\n",
                                                                       frame_index, rgyr_molav[frame_index]);
                     }
/*** If we have used this molecule we have dealt with the first of frame that has to be processed ***/
                     first_of_frame=FALSE;
		   }
		else
		   {
	/*****************************************************/
	/*** Even if not using need to keep track of index ***/
	/*****************************************************/
		      for (this_molecule=0; this_molecule < num_this_mol[itype]; this_molecule++)
											     index++;
		   }
                  start_type=FALSE;
                  printf("End of itype loop\n");
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
	/*********************************************/
	/** Write out in hin file format *************/
	/*********************************************/

			    if (anal_flags[ianal].hin)
			      {
                              }


              printf("End of ianal loop\n");
	    }

          printf("End of frame\n");
	  if (have_config) break;
	 }

/*****************************************************************/
/** Processes requiring information across frames below here *****/
/*****************************************************************/

        if (written_coords) fclose(car_fp);
        if (written_windows) fclose(windows_fp);

        printf("Over entire run sampled %d angles\n", num_angles);

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
/*************************************************************/
/*** Sort out averages for g(r) ******************************/
/*************************************************************/
             /*getting total h-bond lengths and calculating averages for this HISTORY file*/
   if (num_monit >= 0)
      {
         for (iii=0; iii<=num_monit; iii++)
            {
               if ( num_accepted[iii] > 0 ) mean_dist[iii] = cumulat_dist[iii]/num_accepted[iii];
               
               printf("Mean %s %s interaction distance = %10.6f ", monitor_set[iii].atom1,  
                                                                    monitor_set[iii].atom2, 
                                                                    mean_dist[iii] );

               printf("from %d interactions accepted in range %10.6f...%10.6f\n", num_accepted[iii],
                                                                                   monitor_set[iii].min_limit,
                                                                                   monitor_set[iii].max_limit );
/*** This will be moved to after the reading of all frames ***/
/*** Added normalisation factor for g(r) type reporting    ***/

              printf("Bin results for bonds:\n");
              printf("num_accepted[%d] = %d\n", iii,  num_accepted[iii]);
              printf("num_sampled[%d] = %d\n", iii,  num_sampled[iii]);

              for (icheck=0; icheck < monitor_set[iii].bins; icheck++) 
                {
                   distance =  monitor_set[iii].min_limit+icheck*bin_width[iii];

                   distrib[iii][icheck]= cell_volume * distrib[iii][icheck] 
                                                 / (four_pi_delr[iii] *distance *distance*  num_accepted[iii] ) ;

               }
      
/**** Write bond lengths out as csv file ***/
                  sprintf(fname_out, "%s_%s.csv",monitor_set[iii].atom1, monitor_set[iii].atom2);

                  if ( (msd_fp = fopen (&fname_out[0],"w") ) == NULL)
                    {
                      printf ("ERROR Opening %s FILE\n", fname_out);
                      return 1;
                    }

                  write_csv(msd_fp, "distance", "num", "dummy",
                            &bin_distances[iii][0] , &distrib[iii][0], &distrib[iii][0], FALSE,
                             monitor_set[iii].bins-1);

                  printf("Written distance histogram file csv for atoms %s %s, closing file pointer\n",
                                 monitor_set[iii].atom1, monitor_set[iii].atom2);

                  fclose(msd_fp);
             }
	   }

/*** Process angles ***/
   if (num_monit_angle >= 0)
      {
         for (iii=0; iii<=num_monit_angle; iii++)
            {
              if ( num_accepted_angles[iii] > 0)
                   for (icheck=0; icheck < monitor_set_angle[iii].bins; icheck++) 
                        distrib_angle[iii][icheck]= distrib_angle[iii][icheck] /  num_accepted_angles[iii] ;
	      else
                   for (icheck=0; icheck < monitor_set_angle[iii].bins; icheck++) distrib_angle[iii][icheck]= 0.0;
      
/**** Write out angles as csv file ****/
         
                  sprintf(fname_out, "%s_%s_%s.csv",monitor_set_angle[iii].atom1, monitor_set_angle[iii].atom2, monitor_set_angle[iii].atom3);

                  if ( (msd_fp = fopen (&fname_out[0],"w") ) == NULL)
                    {
                      printf ("ERROR Opening %s FILE\n", fname_out);
                      return 1;
                    }

                  write_csv(msd_fp, "Angle", "num", "dummy",
                            &bin_angle[iii][0] , &distrib_angle[iii][0], &distrib_angle[iii][0], FALSE,
                             monitor_set_angle[iii].bins-1);

                  printf("Written Angles csv file closing file pointer\n");

                  fclose(msd_fp);
             }
      }
/*** Process dihedrals ***/
   if (num_monit_dihedral >= 0)
      {
         for (iii=0; iii<=num_monit_dihedral; iii++)
            {
              if ( num_accepted_dihedrals[iii] > 0)
                   for (icheck=0; icheck < monitor_set_dihedral[iii].bins; icheck++) 
                        distrib_dihedral[iii][icheck]= distrib_dihedral[iii][icheck] /  num_accepted_dihedrals[iii] ;
	      else
                   for (icheck=0; icheck < monitor_set_dihedral[iii].bins; icheck++) distrib_dihedral[iii][icheck]= 0.0;
      
/**** Write out dihedrals as csv file ****/
         
                  sprintf(fname_out, "%s_%s_%s_%s.csv", monitor_set_dihedral[iii].atom1,
                                                        monitor_set_dihedral[iii].atom2, 
                                                        monitor_set_dihedral[iii].atom3, 
                                                        monitor_set_dihedral[iii].atom4);

                  if ( (msd_fp = fopen (&fname_out[0],"w") ) == NULL)
                    {
                      printf ("ERROR Opening %s FILE\n", fname_out);
                      return 1;
                    }

                  write_csv(msd_fp, "dihedral", "num", "dummy",
                            &bin_dihedral[iii][0] , &distrib_dihedral[iii][0], &distrib_dihedral[iii][0], FALSE,
                             monitor_set_dihedral[iii].bins-1);

                  printf("Written Dihedrals csv file closing file pointer\n");

                  fclose(msd_fp);
             }
      }

/*** End of csv file writing ***/
/**** End of H-Bond reporting section *****/

             /*calculating averages over all frames*/

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
/*** Output csv file for radius of gyration *****/
/************************************************/
             if (anal_flags[ianal].rgyr)
               {
                  sprintf(fname_out, "rgyr_molav.csv");

                  if ( (msd_fp = fopen (&fname_out[0],"w") ) == NULL)
                    {
                      printf ("ERROR Opening %s FILE for radius of gyration output\n", fname_out);
                      return 1;
                    }

                  write_csv(msd_fp, "time", "rgyr", "dummy",
                            &time_now[0] , &rgyr_molav[0], &rgyr_molav[0], FALSE,
                            frame_index-1);

                  printf("Written rgyr_molav file closing file pointer\n");

                  fclose(msd_fp);

/**** At this point time_now is an array containing the time values for plotting rgyr_molav ***/
/**** rgyr_molav is the radius of gyration listed by frame_index                            ***/
/**** Glib starts                                                                           ***/


/**** Glib ends                                                                             ***/

               }


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
